# Catálogo de nanotubos — esquema do banco e do site

Estado: rascunho para discussão.  A varredura que gera o catálogo já roda no
cluster (vetor `ntb_cat`, fila `normal`, 200 blocos retomáveis).

## 1. O que o catálogo é

Para cada monocamada dos três bancos de origem, todos os tubos distintos com

* diâmetro **t ≤ D ≤ max(30 Å, t + 30 Å)**, onde t é a espessura da camada,
* célula **≤ 50 000 átomos**,
* **a menor célula com resíduo ε ≤ 0,5 %** (a exata quando ela é a menor),

e, para cada sentido de enrolamento, o resultado da checagem de ligações
espúrias da própria ferramenta.

### Por que o diâmetro depende da espessura

Um teto fixo de 30 Å parece razoável até olhar a espessura das camadas: a
mediana é 4,5 Å no C2DB, mas o percentil 90 do 2DMatPedia é **23 Å** e o mais
grosso, um filme de CdTe (`2dm-2644`), tem **55,6 Å**.  Enrolado num tubo de
30 Å, esse filme teria a parede atravessando o eixo.  Medido na primeira
varredura, com o teto fixo:

| | tubos |
|---|---|
| D < t (parede atravessa o eixo) | 586 141 (7,4 %) |
| cavidade interna < 4 Å | 1 168 231 (14,8 %) |
| cavidade interna < 10 Å | 2 471 011 (31,3 %) |

e 1 526 materiais não tinham **nenhum** tubo com cavidade de 10 Å.  Do outro
lado, a fração de tubos sem ligação espúria cresce monotonicamente com o
diâmetro — 28,5 % em D < 10 Å, 73,1 % em 25–30 Å — ou seja, o teto fixo é
justamente onde moram as ligações espúrias: os 1 885 materiais sem nenhum tubo
limpo têm espessura média de 8,8 Å contra 6,1 Å do banco todo.

Daí a regra: piso em D = t, porque abaixo dele o tubo se auto-intersecta, e
teto acompanhando a espessura, para todo sistema poder ter até 30 Å de
cavidade.  A cavidade `d_inner = D − t` fica gravada em cada linha; quem
quiser um canal de 10 Å filtra no site.  Custo medido numa amostra de 496
materiais rodando a regra de verdade: 7,9 milhões de tubos, 15,3 G átomos,
**25 min** em 200 tarefas.

Não guardamos arquivos de estrutura: 7,2 milhões de tubos dariam centenas de
gigabytes, e o construtor monta qualquer um deles em **3 ms**.  O catálogo
guarda os números; o arquivo nasce no pedido.

## 2. Origens e citação

| Origem | Estruturas | Licença / exigência |
|---|---|---|
| C2DB (DTU) | 17 001 | CC BY-NC; citar Haastrup 2018 e Gjerding 2021 |
| 2DMatPedia | 6 351 | citar Zhou et al., Sci. Data 6, 86 (2019) |
| Redes de carbono 2D (jz1c03193) | 1 234 | material suplementar aberto; citar o artigo |

Todo download leva um `CITATION.txt` com as citações das origens presentes na
seleção, e um `manifest.csv` com origem, identificador, (n,m), sentido, ε e a
regra de célula — o suficiente para reproduzir o conjunto.

## 3. Esquema SQLite

```sql
CREATE TABLE materials (
  id INTEGER PRIMARY KEY, source TEXT, uid TEXT, formula TEXT,
  lattice_type TEXT, sector_deg REAL, a REAL, b REAL, gamma REAL,
  n_basis INTEGER, thickness REAL, elements TEXT,          -- " C S Mo "
  citation_key TEXT, UNIQUE(source, uid));

CREATE TABLE tubes (
  id INTEGER PRIMARY KEY, material_id INTEGER REFERENCES materials(id),
  n INTEGER, m INTEGER, diameter REAL, t_norm REAL, atoms INTEGER,
  eps REAL, exact INTEGER, t1 INTEGER, t2 INTEGER, theta REAL,
  d_inner REAL,                                             -- D - espessura
  spurious_ccw TEXT, spurious_cw TEXT,                       -- '' ou 'Se-Se,S-Se'
  clean_ccw INTEGER, clean_cw INTEGER,                       -- 0/1, indexável
  source TEXT, lattice_type TEXT,                            -- repetidos aqui
  UNIQUE(material_id, n, m));

CREATE INDEX tubes_q   ON tubes(clean_ccw, lattice_type, exact, diameter, eps, atoms, d_inner);
CREATE INDEX tubes_src ON tubes(source, clean_ccw, diameter, eps, atoms, d_inner);
CREATE INDEX tubes_all ON tubes(diameter, atoms, eps, exact, clean_ccw, d_inner);
CREATE INDEX tubes_cav ON tubes(d_inner, clean_ccw, atoms, eps);
CREATE INDEX tubes_material ON tubes(material_id, clean_ccw, diameter, eps, atoms);
CREATE INDEX materials_class ON materials(lattice_type, sector_deg);
CREATE TABLE material_elements (material_id INTEGER, element TEXT);
CREATE INDEX melem ON material_elements(element, material_id);
```

Uma linha por (material, n, m), não por sentido de enrolamento: a geometria é a
mesma nos dois sentidos, só as ligações espúrias mudam, e essas vão em duas
colunas.  `material_elements` existe para a busca por elemento ser indexada
("tem Mo", "não tem metal de transição") em vez de varrer texto.

### O que as consultas exigiram (medido no banco de 7,9 M linhas)

O esquema acima tem três coisas que parecem redundância e não são.  Sem elas o
filtro do site não responde:

| consulta | ingênua | com o esquema |
|---|---|---|
| hexagonal, D 5–30 Å, ≤ 5000 átomos, limpo | 18 s (JOIN com `materials`) | **213 ms** |
| contém Mo e é limpo | 12,5 s | **15 ms** |
| sem metais de transição | 42 s (`NOT EXISTS`) | **227 ms** |

1. `source` e `lattice_type` repetidos em `tubes`: sem isso toda contagem vira
   um JOIN de milhões de linhas.
2. `clean_ccw`/`clean_cw` como 0/1: `spurious = ''` não entra em índice.
3. Índice por material **com as colunas do filtro embutidas**: é por aí que
   passa a busca por elemento.
4. Filtro de elemento sempre pela lista de materiais
   (`material_id IN (SELECT id FROM materials WHERE id NOT IN …)`), nunca
   `NOT EXISTS` por tubo — 227 ms contra 42 s.

### O cubo de agregados

A contagem exata pelo índice custa centenas de milissegundos; enquanto o
usuário arrasta um controle, isso pesa.  A tabela `summary` agrega os tubos por
(origem, classe, exato, limpo, D inteiro, faixa de átomos, faixa de ε, faixa de
cavidade) — 20 mil linhas — e responde em **2 ms**, com o mesmo número da
consulta exata (conferido: 1 096 804 pelos dois caminhos).  A estimativa da
caixa de filtro sai do cubo; o número final, do índice.

## 4. Pedidos

```sql
CREATE TABLE requests (
  protocol TEXT PRIMARY KEY,          -- NTB-2026-0912-0042
  email TEXT, name TEXT, filter_json TEXT, filter_hash TEXT,
  n_tubes INTEGER, n_atoms INTEGER, formats TEXT,
  state TEXT,                         -- queued|running|done|failed|expired
  slurm_id TEXT, progress INTEGER, bytes INTEGER,
  created TEXT, finished TEXT, expires TEXT, path TEXT);
CREATE INDEX requests_hash ON requests(filter_hash, state);
```

## 5. Fluxo do pedido

1. **Filtro.**  A caixa de filtro consulta o SQLite e responde em
   milissegundos: quantos tubos, quantos átomos, tamanho estimado (medido:
   35 kB por tubo em XYZ, 11 kB compactado) e tempo estimado (3 ms por tubo).
2. **Confirmação** com e-mail e nome.  Sem teto de tamanho.
3. **Protocolo** gerado na hora, mostrado na tela e enviado por e-mail.
4. **Execução no cluster, nunca no headnode.**  O site grava o pedido e
   submete um `sbatch` na fila **`high`** (os pedidos rodam em minutos).  O
   trabalho constrói os tubos, compacta e devolve para `~/projects`.
5. **Acompanhamento** pelo protocolo: a página lê `progress` e mostra a barra.
   O trabalho atualiza a linha do pedido a cada N tubos.
6. **E-mail com o link** ao terminar, válido por **7 dias**; limpeza
   automática depois disso.
7. **Cache por `filter_hash`**: pedido idêntico ainda válido devolve o mesmo
   arquivo na hora.

## 6. O que o site entrega

* **Um tubo:** construção na hora, em XYZ, CIF, POSCAR ou LAMMPS.
* **Uma seleção:** ZIP com os arquivos, `manifest.csv` e `CITATION.txt`.
* **O catálogo:** o SQLite para baixar, mais um script de dez linhas que
  reconstrói qualquer subconjunto na máquina do usuário.
* **Um pacote curado:** os tubos exatos pequenos (até 200 átomos e 20 Å), que
  é o que a maioria vai querer rodar em DFT.

## 7. Reprodutibilidade

Cada linha do catálogo é reconstruível a partir de (origem, identificador,
n, m, sentido) com a versão da ferramenta e os parâmetros da regra, que ficam
gravados numa tabela `meta` do próprio SQLite.  Nada no catálogo depende de um
limite de busca arbitrário: a célula é a menor dentro de ε ≤ 0,5 %.

## 8. O site: páginas, API e execução

O site já existente (FastAPI em `web/api/main.py`, página única em
`web/static/index.html`) ganha uma aba **Catálogo** ao lado das que já
existem.  Nada do que está lá muda: o catálogo é um módulo novo
(`web/api/catalogue.py`) montado em `/api/cat/*`, com o seu próprio SQLite em
`web/data/catalogue.db` — só leitura para o site, escrito pela varredura.

### 8.1 A aba Catálogo, de cima para baixo

```
┌─ Catálogo de nanotubos ────────────────────────────────────────────────┐
│  24 586 monocamadas · 7,2 M tubos · D ≤ 30 Å · ε ≤ 0,5 %              │
│                                                                        │
│  Origem     [x] C2DB  [x] 2DMatPedia  [x] Redes de carbono            │
│  Elementos  contém [Mo][S]        sem  [  ]      espécies [1..9]      │
│  Rede       [x] hexagonal [x] quadrada [x] retangular [x] oblíqua     │
│  Diâmetro   [ 5 ] a [ 30 ] Å        Átomos/célula  ≤ [ 50000 ]        │
│  Resíduo    ( ) só exatos  (•) ε ≤ [0,5] %                            │
│  Quiralidade [x] armchair [x] zigzag [x] quiral                       │
│  Ligações   (•) sem espúrias  ( ) só com  ( ) todas   sentido [ambos]  │
│  Formato    [x] XYZ [ ] CIF [ ] POSCAR [ ] LAMMPS    vácuo [10] Å     │
│                                                                        │
│  ▸ 12 480 tubos · 3,1 M átomos · ~140 MB compactado · ~2 min          │
│                             [ Ver amostra ]  [ Pedir este conjunto ]  │
└────────────────────────────────────────────────────────────────────────┘
```

Cada mexida num controle refaz a contagem (`POST /api/cat/count`, um
`SELECT COUNT(*)`/`SUM(atoms)` sobre os índices; medido em milissegundos).
"Ver amostra" mostra 12 tubos da seleção, com D, ε, átomos e o mapa de
quiralidade do material — e cada um é construído na hora pelo endpoint de tubo
único, que é o que o site já faz hoje.

O tamanho e o tempo aparecem **antes** do pedido, e o botão de pedir abre a
confirmação com e-mail e nome.  Sem teto: se a seleção der 8 GB, o site diz
"8 GB, ~40 min" e deixa o usuário decidir.

### 8.2 Endpoints

| Método | Rota | O que faz |
|---|---|---|
| `GET`  | `/api/cat/meta` | totais, versão da ferramenta, parâmetros da regra |
| `POST` | `/api/cat/count` | filtro → (tubos, átomos, bytes, segundos) |
| `POST` | `/api/cat/sample` | filtro → 12 linhas, para a amostra |
| `POST` | `/api/cat/request` | filtro + e-mail → protocolo, submete o `sbatch` |
| `GET`  | `/api/cat/status/{protocolo}` | estado, progresso, link quando pronto |
| `GET`  | `/api/cat/download/{protocolo}` | o ZIP, enquanto não expirar |
| `GET`  | `/api/cat/tube/...` | um tubo, construído na hora (já existe) |

`/api/cat/count` e `/api/cat/sample` são consultas puras — o site responde
sozinho.  Só `/api/cat/request` toca o cluster.

### 8.3 Execução, nunca no headnode

`POST /api/cat/request` faz três coisas e devolve: grava a linha em
`requests`, escreve `filter.json` no diretório do pedido e chama

```
sbatch -p high -J ntb_req_<protocolo> -n 1 -c 4 -t 01:00:00 \
       --export=PROTO=<protocolo> <repo>/web/deploy/request.sbatch
```

`request.sbatch` carrega `module load python/3.12`, constrói os tubos do
filtro num diretório do nó, compacta, copia o ZIP para
`~/projects/ntbuilder_cat/out/<protocolo>.zip` e marca `state='done'`.  A cada
200 tubos ele escreve `progress` na tabela — a barra da página é esse número.
Pedido grande é fatiado em vetor (`-a 0-N`), um job de junção depois; a barra
soma os blocos.

O headnode só faz consulta SQL e `sbatch`.  Nenhuma construção roda nele.

### 8.4 E-mail

Dois e-mails por pedido, do `lab.nanoeng@gmail.com`:

1. **na hora**: o protocolo, o filtro em texto, o link de acompanhamento;
2. **ao terminar**: o link do ZIP e a data de expiração (7 dias).

O envio sai de uma fila em disco (`web/data/mail/`), não do request handler:
se o Gmail recusar, o pedido não falha, e o reenvio é uma linha de comando.
As credenciais (senha de app) ficam **fora do repositório**, em
`/etc/ntbuilder/mail.env` lido pelo systemd — nunca em `git`.

### 8.5 Limpeza e limites

* `cleanup_tmp.sh` ganha a varredura dos ZIPs: apagar o que passou de 7 dias
  e marcar `state='expired'`.
* Um pedido idêntico (`filter_hash`) ainda válido devolve o mesmo arquivo sem
  submeter nada.
* Um e-mail com mais de 3 pedidos em fila espera o primeiro terminar — não é
  teto de tamanho, é fila justa.

## 9. Correções de 12/09/2026, antes da versão pública

### Unidade de contagem e sentido de enrolamento

Uma linha de `tubes` é um (n,m) de um Sistema 2D. Ela vale um tubo quando os dois sentidos de enrolamento dão o mesmo tubo e dois quando dão tubos diferentes. O teste é de simetria: enrolar para dentro é enrolar para fora a camada virada (z → −z), e os dois sentidos coincidem quando alguma operação no plano leva a camada virada na original e (n,m) na sua própria órbita. Medido nos exemplos, grafeno, MoS₂ e bifenileno têm sempre um sentido, o MoSSe tem sempre dois (50 de 50) e o penta-grafeno tem dois em 82 de 92 (n,m). As colunas `senses`, `n_all` e `n_clean` guardam isso, e toda contagem soma `n_all` ou `n_clean`, nunca conta linhas. O pedido entrega os dois arquivos quando os sentidos diferem, e o usuário não escolhe sentido.

A tolerância desse teste é geométrica, 0,02 Å no plano e na altura, e não a da busca de simetria. Com 1e-3 Å na altura, 7 % dos Sistemas 2D recebiam dois sentidos sem precisar.

### Checagem de ligações espúrias periódica

`check_spurious_bonds` procurava pares dentro de uma célula só do tubo. Um par próximo que atravessa a borda da célula passava despercebido, e a borda cai em lugar diferente conforme o corte. O catálogo revelou isso: 2 727 tubos de sentido único tinham vereditos diferentes nos dois sentidos, e 58 de 60 deles concordam quando as células vizinhas entram na busca. A checagem agora inclui imagens ao longo do eixo até cobrir o maior corte. Numa amostra de 400 tubos nenhum veredito mudou, e nos 276 pontos do `mosse_senses` do artigo também não.

### Tolerância de simetria

A tolerância das posições atômicas era 5 vezes a tolerância relativa da rede, em fração da célula, e chegava a 1,18 Å (2,2 % do catálogo acima de 0,5 Å). Passou a ser uma escada em Å, 0,001, 0,01, 0,1 e 0,5 Å, como no Materials Studio, com o degrau mais apertado que dá o maior grupo pontual. A escada da rede continua relativa (até 2 %), porque aplicar 0,5 Å aos parâmetros de rede transformou o bifenileno em quadrado e promoveu 321 de 3 000 sistemas. Numa amostra de 3 000, nenhuma classe muda e 64 setores mudam; os Sistemas 2D das figuras aprovadas fecham a 0,001 ou 0,01 Å e não mudam.

### Quiralidade

`kind` é 0 para aquiral e 1 para quiral, pelo grupo pontual: aquiral quando uma reflexão U do grupo leva (n,m) em ±(n,m). A classificação anterior por índices chamava de armchair o (n,n) do bifenileno e de zigzag o (n,0) do penta-grafeno, e os dois são quirais.

### Filtro de espessura e cubo

O cubo ganhou a dimensão `t_bin` (bordas 0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30 e 60 Å) e perdeu `clean_ccw`/`clean_cw`, trocados pelas somas `n_clean` e `n_atoms_clean`.

### Defeito conhecido, ainda aberto

`core/io.py::lattice_type_at` chama de oblíqua a rede retangular centrada quando a célula reduzida tem a ≠ b (2b·cosγ = ±a). No catálogo são 1 167 Sistemas 2D, e 372 têm espelho de fato. A quiralidade e o sentido não dependem do rótulo, só o filtro de classe de rede.

## 10. Correções e bases novas de 13/09/2026

### Critério de ligação espúria com encurtamento mínimo

Um par de espécies conta como ligação espúria quando, no tubo, fica abaixo do corte de 1,2 vez a soma dos raios covalentes sem ser ligado na camada plana, **e** quando a sua menor distância no tubo é pelo menos 10 % menor que a menor distância do mesmo par na camada plana. A segunda condição é uma histerese sobre o corte. Sem ela, um par que já nasce logo acima do corte é marcado por qualquer curvatura. O caso que revelou isso é o Fe₂Mo₂F₂O₈ do C2DB (`2FFeMoO4-1`), onde o par Fe–Mo atravessa a ponte Fe–O–Mo a 3,71 Å com corte de 3,67 Å, e todos os 125 tubos eram marcados, inclusive o (17,2) de 36,6 Å por um aperto de 3,5 %. Com 10 %, 115 desses 125 tubos ficam limpos e os 10 que continuam marcados têm D ≤ 13 Å. Numa amostra de 250 tubos do catálogo com alerta, a mediana de encurtamento dos pares marcados era 26 %, e 10,4 % dos tubos perdem o alerta. Os vereditos das figuras do artigo não mudam: os 552 do `mosse_senses` e os 836 do `exotic_spurious` (418 pontos nos dois sentidos, células de até 2,7 milhões de átomos, conferidos no cluster). O parâmetro é `min_shortening` em `check_spurious_bonds`, e zero devolve o teste só pelo corte. Na mesma noite o catálogo passou para a checagem par a par da seção 14, que marca também as ligações rompidas.

### Células fora do plano

1 628 dos 6 351 Sistemas 2D do 2DMatPedia chegam com a ou b fora do plano xy, às vezes com o vetor de vácuo entre eles (até 22,6 Å em z). A conversão original usava só x e y desses vetores e distorcia a rede, com mudança de área de uns 2 % e de espessura de até 1 Å. A conversão agora gira célula e átomos juntos até o plano da camada ficar em xy, o que preserva todas as distâncias (diferença de 10⁻¹⁵ Å no teste). O C2DB não tem nenhum caso, e o JARVIS tem 160.

### Repetidos entre bases

Repetido é quem tem a mesma fórmula reduzida, área por átomo e espessura próximas (pré-filtro), o mesmo número de átomos na célula primitiva e, decisivo, as duas estruturas sobrepostas com todo átomo perto do vizinho de mesma espécie. A sobreposição procura a correspondência entre as redes (matriz inteira M com |det M| igual à razão entre os números de átomos, o que aceita supercélulas e outra escolha de base), tira a escala isotrópica da rede, tenta a camada virada em z e cada translação, e mede o maior deslocamento de um átomo em unidades de √(área por átomo). Dentro da mesma base o limite é 0,10, e entre bases diferentes é 0,20 (cerca de 0,3 Å), porque ali é o mesmo composto relaxado por outro grupo. Quem entra como apelido de um Sistema 2D ainda tem de ser igual, pela regra da mesma base, a todo apelido da sua base que já está nele.

A primeira versão usava uma impressão digital de distâncias com corte de 2,6 vezes a escala. Ela separava o MoS₂ 2H do 1T, mas não via o empilhamento, porque entre camadas afastadas a distância 3D passa do corte e é justamente ali que ABA e ABC diferem. Resultado medido, 609 geometrias distintas da mesma base juntadas, quase todas MXenes. A calibração da sobreposição (`dedup_cal/cal2.py`) deu estes números. Nos 10 956 pares de geometrias da mesma família de MXenes o menor desvio é 0,97, e nenhum fica abaixo de 0,3, nem os de mesma energia, de modo que as 4 489 MXenes são todas estruturas distintas. Nos apelidos da mesma base da versão anterior, 56 ficam perto de zero e 609 acima de 1,6, sem nada entre 0,05 e 0,3. Nas 295 entradas do MatHub-2d declaradas cópia do JARVIS, 267 ficam até 0,2, contra 224 reconhecidas pela impressão digital. Nos apelidos entre bases da versão anterior a mediana é 0,004 e o percentil 95 é 0,08, e 184 pares acima de 0,15 eram estruturas diferentes.

### Bases novas

| Base | Lidas | Repetidas | Novas | Licença |
|---|---|---|---|---|
| MC2D (Mounet 2018, Campi 2023) | 2 742 | 990 | 1 752 | CC BY 4.0 |
| JARVIS-DFT 2D (Choudhary 2017, 2020) | 1 103 | 817 | 286 | CC BY 4.0 |
| MXenes (Ontiveros 2025) | 4 489 | 301 | 4 188 | MIT |
| Alexandria 2D até 0,1 eV/átomo acima do hull (Wang 2023, Cavignac 2026) | 19 691 | 2 954 | 16 737 | CC BY 4.0 |
| MatHub-2d (Yao 2023) | 1 907 | 1 661 | 246 | não declarada; 4 estruturas sem rede relaxada ficaram de fora |

Das 4 489 MXenes nenhuma repete outra MXene, e as 301 repetidas estão no C2DB (299) e no MC2D (2). O MatHub-2d declara 1 517 estruturas tiradas do C2DB, e o critério geométrico reconhece 1 518 repetidas do C2DB, uma verificação que não entrou na calibração. A energia acima do hull do Alexandria e a energia relativa das MXenes ficam gravadas.

### Catálogo final (13/09/2026, camadas corrigidas e um tubo por órbita)

46 403 Sistemas 2D principais e 8 093 apelidos, 14 609 330 pares (n,m), 20 334 372 nanotubos (sentidos distintos contam dois), 11 608 889 limpos (57.1 %) pela checagem par a par da seção 14. Pela checagem anterior, só de ligações novas, eram 14 927 476 limpos (73.4 %). Nenhum Sistema 2D ficou sem tubo, nenhuma camada guardada tem vão interno maior que 6 Å, nenhum Sistema 2D de uma amostra de 300 tem dois tubos na mesma órbita de simetria, e o total do cubo bate com a soma direta da tabela `tubes`.

| Base | Sistemas 2D | Nanotubos | Limpos |
|---|---|---|---|
| C2DB | 16 990 | 9 638 917 | 5 607 433 |
| Alexandria 2D | 16 713 | 5 037 658 | 3 629 906 |
| 2DMatPedia | 5 079 | 2 484 388 | 1 211 562 |
| MXenes | 4 172 | 2 169 374 | 588 015 |
| MC2D | 1 744 | 683 071 | 338 150 |
| Redes de carbono 2D | 1 232 | 119 451 | 119 132 |
| MatHub-2d | 243 | 107 229 | 55 328 |
| JARVIS-DFT 2D | 230 | 94 284 | 59 363 |

Antes da correção das camadas partidas (seção 12) eram 46 905 Sistemas 2D e 22 707 013 nanotubos, e antes de agrupar os (n,m) equivalentes (seção 13) eram 21 859 314 nanotubos. A primeira diferença vem dos 1 323 Sistemas 2D varridos de novo com a camada inteira (617 698 linhas antigas descartadas na carga), das 22 camadas recusadas por defeito de origem e das cópias do C2DB que o dedup passou a reconhecer. A segunda vem de 1 011 060 pares (n,m) que eram cópias por simetria de outros do mesmo Sistema 2D.

Sobreposição entre bases (Sistemas 2D presentes nas duas) em `overlap.json`. 42 408 Sistemas 2D aparecem numa só base, 3 306 em duas e 1 191 em três ou mais.

Nenhum (Sistema 2D, n, m) aparece duas vezes nos shards. O `load_db.py` descarta repetidos por segurança e contou zero. Uma contagem paralela (`count_nts.py`) chegou a dar 187 516 pares a menos porque guardava só `hash((uid, n, m))`, e no CPython `hash(-1) == hash(-2)`, o que juntava (n,−1) e (n,−2) nas redes com índices negativos. Conferido nas linhas brutas do C2DB `2Ge-2` (1 345 pares distintos, os mesmos do banco), e o script agora guarda a chave inteira. Tubos por Sistema 2D caem com a área da célula, de 430 a 530 abaixo de 30 Å² para 67 acima de 120 Å², e sobem com a espessura, porque a janela t ≤ D ≤ t + 30 Å cresce com t.

## 11. Página e entrega (13/09/2026)

### Entrega na hora ou por e-mail

A construção é estimada por t = 1,1 × 10⁻⁵ s por átomo + 1,3 × 10⁻³ s por tubo, vezes o número de formatos. Os coeficientes vêm de medida no headnode, com um núcleo, nice 10 e CIF: 1,2 s para 257 tubos pequenos do C2DB (82 mil átomos), 3,8 s para 300 MXenes (351 mil) e 4,4 s para 209 tubos da seleção padrão (388 mil). O ajuste reproduz os três com erro abaixo de 10 %, e um pedido de teste com estimativa de 19,6 s levou 19 s. Os 3 ms por tubo usados antes subestimavam de 2 a 7 vezes.

Até 60 s estimados (`NTB_INSTANT_MAX_SECONDS`), o pedido é montado no próprio servidor com o mesmo `request_job.py` do cluster (a cópia do NFS, `site/`), que escreve `progress.json` e `done.json` do mesmo jeito. A página mostra uma barra de progresso e dispara o download ao terminar, e o observador de pedidos manda o e-mail com o protocolo e o link de 7 dias. Acima disso, ou com mais de `NTB_INSTANT_QUEUE_MAX` (6) montagens esperando, o pedido vai para a fila `high` e o link chega por e-mail, como antes.

Reserva do headnode: `nice -n 10`, teto de memória de `NTB_INSTANT_MEM_GB` (6 GB), `NTB_INSTANT_SLOTS` (1) montagem por worker do uvicorn (2 workers, então no máximo 2 ao mesmo tempo) e fila de espera limitada. Uma montagem local que morre sem escrever o resultado (reinício do serviço, falta de memória) vira falha com e-mail, e uma que nunca começou falha depois de 1 h. O estado passa para pronto ou falha com `UPDATE ... WHERE state IN ('queued','running')`, e só quem de fato mudou a linha manda o e-mail, porque com dois workers a consulta do protocolo e a ronda chegavam a disputar o mesmo pedido.

O protocolo do dia avança até achar uma pasta de saída que não existe, e a pasta é criada com `exist_ok=False`. Um servidor de teste com outro banco de pedidos numerou um pedido como `NTB-20260913-0001`, a mesma pasta de um pedido real no NFS, e sobrescreveu o arquivo dele. O pedido foi refeito idêntico (9 004 674 bytes) a partir do `filter_json` guardado e do banco que estava no ar na hora.

### Amostra paginada

A amostra fica sempre aberta e acompanha o filtro. Cada linha é um (n,m) de um Sistema 2D, com os arquivos dos sentidos que passam no filtro, e com "Sem ligação espúria" a coluna de espúrias some. A ordem é a do banco: `ORDER BY atoms` custava 2,5 s na seleção padrão (árvore temporária sobre 11,6 milhões de linhas), e `ORDER BY id` responde em milissegundos em qualquer página (0,6 s na página 100 001). O número de páginas vem de `n_rows` e `n_rows_clean`, duas colunas novas do cubo `summary` (conferidas contra a tabela, 16 005 587 e 11 569 559), e com filtro de elementos vem do `COUNT(*)` que acompanha a soma no índice.

### Classes de elementos

"Incluir classes" pede ao menos um elemento de cada classe marcada (`elements_any`, uma lista de grupos, cada grupo vira `id IN (SELECT material_id FROM material_elements WHERE element IN (...))`), e "Excluir classes" marca os elementos dela como excluídos. O cubo não responde filtro de elementos, então a contagem sai do índice (2,5 s para "com metal de transição" sobre 7,8 milhões de linhas).

### Layout

Topo com a apresentação e a tabela periódica em 4/5 e a Base de Dados 2D em 1/5, filtro em 1/3 ao lado da contagem e da amostra em 2/3, pedido e explicação na largura toda, tudo em uma coluna abaixo de 1000 px. Os balões da explicação ficam em colunas de larguras diferentes. A página mede a altura real de cada balão em todas as larguras, de 4 em 4 px, numa passada só (todos os clones fora da tela e um único cálculo de layout, contra 3 s lendo altura por altura), testa as 99 divisões dos 8 balões em até 5 colunas na ordem de leitura, acha para cada uma a menor altura em que as colunas cabem na largura, e fica com a menor soma de altura e sobra média por coluna. Medido no Chromium sem tela, a sobra entre colunas fica entre 0 e 20 px (uma linha de texto repartida) em janelas de 400 a 1440 px, e o cálculo leva de 0,1 a 0,6 s. Uma grade densa com balões de uma ou duas colunas, tentada antes, deixava três linhas com espaço vazio. A API anuncia `features` no `/meta`, e a página esconde classes e páginas quando falam com uma API que não as conhece.

## 12. Camadas partidas pela borda da célula (13/09/2026)

Nenhum conversor juntava os átomos pela periodicidade ao longo do vácuo. As bases guardam as posições dentro da célula, e uma camada centrada perto de c = 0 chegava partida em duas metades separadas por quase uma célula inteira. No 2DMatPedia, o 2dm-1 (IrF₂) tem os átomos em c = 0,0, 0,966 e 0,034, e o manifesto registrava 22 Å de espessura para uma camada de uns 2 Å. Com um vão interno maior que 8 Å havia 1 654 Sistemas 2D no banco (1 518 do 2DMatPedia, 65 do JARVIS, 65 MXenes, 5 do Alexandria, 1 do MatHub-2d), de onde vinham 1,36 milhão de nanotubos construídos a partir de camadas quebradas, com espessura, janela de diâmetro e ligações espúrias erradas. C2DB, MC2D e as redes de carbono não tinham nenhum caso. O defeito apareceu ao abrir no construtor o CIF de um Sistema 2D do catálogo, quando o leitor de CIF separou as duas metades.

`make_manifest.unwrap_slab` acha o eixo de vácuo pelo mesmo critério da rotação para o plano, toma o maior vão periódico das frações ao longo dele como o vácuo e desce uma célula os átomos acima desse vão. A função entra em `to_layer_frame` e nos ramos sem rotação do C2DB e do 2DMatPedia. Depois da correção, o número de camadas com vão interno maior que 6 Å caiu de 1 519 para 8 no 2DMatPedia, de 65 para 1 no JARVIS e de 66 para 0 nas MXenes.

O que sobra são defeitos das próprias bases, e `row()` os recusa (vão interno maior que 6 Å, bem acima da distância de van der Waals). São duas folhas na mesma célula (2dm-2924 AgBi₃, 2dm-2933 GaAs, JVASP-60475 BaB₂Se₆), um átomo solto a 9 Å (2dm-2510 Hf₂O₇, 2dm-2607 Zr₂O₇) e três átomos espalhados pela célula (2dm-6043 LiSn₂, 2dm-6060 RbSn₂, 2dm-6098 NaPb₂, 2dm-6166 KSn₂). Dentro do corte de 0,1 eV/átomo o Alexandria perdeu 13, conferidas uma a uma. Oito são dois planos atômicos planos a 6,1–6,7 Å um do outro, sem ligação entre eles (agm2000036445 Au₄Br₄, agm2000025052 Cr₂Cl₄, agm2000025189 Cu₂Br₄, agm2000022246 Cu₄Cl₄, agm2000019546 Hg₂Cl₄, agm2000022491 Hg₂I₄, agm2000002145 I₄H₄, agm2000019598 Pt₂Cl₄). Três são duas camadas de haleto de lantanídeo empilhadas na mesma célula, com 6,3 a 8,1 Å de vão (agm2000033178 Ho₄Cl₁₂, agm2000022230 Sm₂Cl₆, agm2000033240 Tb₄Cl₁₂). Duas são um plano só de hidrogênio a 6–7 Å dos metais (agm2000001555 H₆Tl₂, agm2000001089 H₆Zn₃).

Com as camadas inteiras, o dedup reconhece mais cópias reais: 1 217 estruturas do 2DMatPedia repetem o C2DB, contra 856 antes. O manifesto final tem 46 403 Sistemas 2D e 8 093 apelidos, sem nenhuma estrutura distinta engolida. `changed_systems.py` compara cada Sistema 2D com a camada que de fato foi varrida e separou 1 323 cuja camada mudou (1 143 do 2DMatPedia, 113 do Alexandria, 50 MXenes, 12 do JARVIS, 4 do MC2D, 1 do C2DB). Só eles foram varridos de novo (`catalogue_fix.sbatch`, shards `shard-W*`), e a carga (`NTB_REPLACED`, `NTB_REPLACED_PREFIX`) descarta as linhas antigas desses sistemas vindas de qualquer outro shard.

Verificação para qualquer base nova: contar os registros cujo maior vão interno em z passa de 6 Å.

## 13. (n,m) equivalentes guardados como tubos distintos (13/09/2026)

`core.planegroup.unique_indices`, que escolhe um (n,m) por família de índices equivalentes, tinha `tol=2e-3` fixo como padrão, enquanto `chirality_group` e `sector_deg` usam a tolerância da busca de simetria (`structure.sym_tol`, que a escada de tolerâncias em Å dos átomos ajusta). Numa camada levemente distorcida as duas decisões discordavam. O C2DB 1AlC-1 tem γ = 119,71°, e depois do ajuste os átomos ficam em (0,349, 0,318) e (0,680, 0,653), perto de (1/3, 1/3) e (2/3, 2/3) para a tolerância do ajuste, mas não para 2e-3. O setor saía com 12 operações (30°) e a redução dos índices com 4 (90°). No construtor, o mapa desenhava o setor de 30° e espalhava 731 pontos por 90°, que o usuário notou na página. No catálogo, a varredura (`scan_chirality(unique_only=True)`) guardou os (n,m) equivalentes como tubos distintos. Numa amostra de 600 Sistemas 2D, 45 (7,5 %) tinham o grupo reduzido pela metade, do C2DB, do 2DMatPedia e do Alexandria.

O padrão de `unique_indices` passou a ser `tol=None`, a mesma tolerância das outras duas funções. Nos 14 exemplos do construtor, que são as camadas das figuras do artigo, as duas tolerâncias dão exatamente o mesmo resultado, e os 113 testes do núcleo passam. O 1AlC-1 fica com 258 pontos, todos dentro do setor de 30°.

Não foi preciso varrer de novo. Os pares a mais já estavam no banco como cópias por simetria. `load_db.py` calcula o grupo de cada Sistema 2D com a tolerância do ajuste, agrupa as linhas de cada Sistema 2D pelas órbitas desse grupo e guarda uma por órbita, com a regra de nome de `unique_indices` (m ≥ 0 e o maior n). Saíram 1 011 060 pares, 6,5 % dos 15 620 390, e o 1AlC-1 ficou com os mesmos 258 (n,m) que o construtor dá até D = 29,9 Å.

## 14. Ligações formadas e rompidas, par a par (13/09/2026, noite)

A checagem de ligações espúrias via só metade do efeito da curvatura. Ela marcava um par de espécies que, no tubo, ficava abaixo do corte sem se ligar em nenhum ponto da camada plana. Não via duas coisas. A primeira é a ligação rompida, quando um par ligado na camada passa do corte no tubo. A segunda é a ligação nova de um par de espécies que já se liga em outro ponto da camada. O grafeno (2,0), de D = 1,57 Å, tem as paredes opostas mais perto que o corte C–C e passava limpo, porque C–C já é ligado no plano.

A checagem nova é `core.builder.check_curvature_bonds`. Cada átomo do tubo guarda de qual átomo da célula 2D e de qual célula da rede ele vem (`_iter_sites`), então todo par de átomos tem uma distância na camada e outra no tubo, com o mesmo corte de `compute_bonds` (1,2 vez a soma dos raios covalentes).

* Formada: abaixo do corte no tubo, acima dele na camada, e ao menos 10 % mais curta que na camada. Um par que no tubo fica mais perto que a distância mínima de ligação também conta.
* Rompida: abaixo do corte na camada, acima dele no tubo, e ao menos 10 % mais longa que na camada.

A margem de 10 % é a histerese da seção 10, agora por par de átomos e nos dois sentidos. As ligações da camada são seguidas até o tubo de forma analítica: o parceiro de deslocamento (dx, dy) na rede está em (s + ds, w + dw) no tubo, com ds = (dx·t2 − dy·t1)/D e dw = (n·dy − m·dx)/D, sem busca de vizinhos na costura. As formadas vêm de uma busca de vizinhos no tubo com as imagens axiais, e cada par volta à camada somando a translação T por imagem e Ch por volta. No banco, `spurious_ccw`/`spurious_cw` guardam `+A-B` para a ligação formada e `-A-B` para a rompida, e a página escreve "formada" e "rompida". `check_spurious_bonds` continua no código para a interface desktop e para as figuras do artigo, cujos vereditos (`mosse_senses`, `exotic_spurious`) foram calculados com ela.

Conferências. As posições que a checagem usa são as de `build_nanotube` (diferença abaixo de 10⁻⁶ Å). Numa amostra de 1 500 linhas sorteadas (2 104 tubos), todo par marcado pela checagem antiga continua marcado, sem exceção. No ClOP do C2DB, tubo (9,6) de D = 13,4 Å, a ligação P–O de 2,008 Å na camada vai a 2,221 Å no tubo (+10,6 %, corte 2,076 Å). No Mo₂CCl₂ ABA_H, tubo (16,16), nenhum átomo perde ou ganha vizinho nos dois sentidos, e os vãos que o visualizador mostra na ponta da célula são só ligações não desenhadas através da borda.

Por que tanta ligação rompida em diâmetro grande. Enrolar sem relaxar estica a face externa na razão (R + z)/R e comprime a interna. Numa camada de 11 Å num tubo de 39 Å a face externa estica 28 %, e a C–H do Rb₂C₆H₆S₆ do C2DB passa do corte. Por isso as MXenes, grossas, ficam com 27,1 % dos tubos limpos, e as redes de carbono 2D, planas, com 99,7 %.

Campanha. `recheck.py` monta cada tubo com o (n, m, t1, t2) gravado, sem refazer a busca, e checa os dois sentidos com os mesmos sítios. Foram 150 tarefas na fila `normal` (`recheck.sbatch`, uns 22 min cada, 21 µs por átomo), lendo uma cópia do banco em modo imutável, sem nenhum erro. `apply_recheck.py` grava numa cópia os resultados das 14 609 330 linhas, pela regra do `load_db` (sentido único só é limpo se os dois sentidos concordam), refaz `summary` e `material_cube` e confere as somas com a tabela. 3 192 linhas de sentido único tiveram resultados diferentes nos dois sentidos e ficaram sem o selo.

Resultado: 20 334 372 nanotubos, 11 608 889 limpos (57.1 %), antes 14 927 476 (73.4 %). Linhas (n,m) com algum tubo limpo: 8 529 169.

| Base | Nanotubos | Limpos | Antes |
|---|---|---|---|
| C2DB | 9 638 917 | 5 607 433 (58.2 %) | 6 937 832 |
| Alexandria 2D | 5 037 658 | 3 629 906 (72.1 %) | 4 130 780 |
| 2DMatPedia | 2 484 388 | 1 211 562 (48.8 %) | 1 880 927 |
| MXenes | 2 169 374 | 588 015 (27.1 %) | 1 191 008 |
| MC2D | 683 071 | 338 150 (49.5 %) | 506 253 |
| Redes de carbono 2D | 119 451 | 119 132 (99.7 %) | 119 451 |
| MatHub-2d | 107 229 | 55 328 (51.6 %) | 81 698 |
| JARVIS-DFT 2D | 94 284 | 59 363 (63.0 %) | 79 527 |
