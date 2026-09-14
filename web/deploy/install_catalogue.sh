#!/bin/bash
# Instala o catálogo para o site e para o cluster, em duas etapas.
#
#   install_catalogue.sh          código do job no NFS + banco novo PREPARADO ao
#                                 lado do que está no ar (catalogue.db.new, com
#                                 ANALYZE), sem tocar no site
#   install_catalogue.sh --swap   troca o banco do site pelo preparado (mv
#                                 atômico) e devolve a cópia com ANALYZE ao NFS
#
#   web/data/catalogue.db   cópia local, que o site lê (rápido, sem NFS)
#   ~/projects/ntbuilder_cat/site/   o que o job do nó precisa: core/, o
#                                    tradutor do filtro, o construtor e o job
#
# Por que duas etapas: o site abre uma conexão por requisição, e um cp por cima
# do arquivo aberto entrega a quem estiver lendo um banco pela metade.  Entre
# as duas etapas dá para testar o banco novo numa porta separada
# (NTB_CATALOGUE_DB=web/data/catalogue.db.new).
#
# O nó não vê /home do headnode (é local de cada máquina), então o código do
# job tem de estar no NFS.  Rode isto depois de cada mudança em core/ ou nos
# módulos do catálogo, senão o pedido roda com código velho.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
NFS=${NTB_PROJECTS:-$HOME/projects}/ntbuilder_cat
SITE=$NFS/site
LIVE=$REPO/web/data/catalogue.db
NEW=$LIVE.new

if [ "${1:-}" = "--swap" ]; then
    [ -f "$NEW" ] || { echo "ERRO: $NEW não existe — rode primeiro sem --swap"; exit 1; }
    [ -f "$LIVE" ] && ln -f "$LIVE" "$LIVE.anterior"     # o banco que sai fica guardado
    mv -f "$NEW" "$LIVE"
    cp "$LIVE" "$NFS/catalogue.db.tmp" && mv -f "$NFS/catalogue.db.tmp" "$NFS/catalogue.db"
    ls -la "$LIVE"
    echo "banco trocado; o anterior ficou em $LIVE.anterior"
    exit 0
fi

mkdir -p "$SITE" "$NFS/out" "$NFS/logs" "$REPO/web/data"
rsync -a --delete --exclude __pycache__ "$REPO/core/" "$SITE/core/"
cp "$REPO/web/api/cat_build.py" "$REPO/web/api/cat_query.py" "$SITE/"
cp "$REPO/web/deploy/request_job.py" "$SITE/"
echo "código do job em $SITE"

if [ -f "$NFS/catalogue.db" ]; then
    cp "$NFS/catalogue.db" "$NEW"
    # ANALYZE e o que faz o planejador escolher o indice certo no filtro por
    # elemento: sem as estatisticas, "contem Mo e S" custava 10 s em vez de
    # 176 ms.  Roda em 33 s.
    sqlite3 "$NEW" "ANALYZE;" 2>/dev/null \
        || python3 -c "import sqlite3,sys; c=sqlite3.connect(sys.argv[1]); c.execute('ANALYZE'); c.commit()" "$NEW"
    ls -la "$NEW"
    echo "banco novo preparado em $NEW; teste e depois rode com --swap"
else
    echo "AVISO: $NFS/catalogue.db não existe — rode a varredura e load_db.py"
fi
