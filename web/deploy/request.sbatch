#!/bin/bash
# Um pedido do catálogo, na fila rápida.  O site submete isto; o headnode só
# consulta SQL e chama sbatch.
#SBATCH -J ntb_req
#SBATCH -p high
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -o /mnt/home/mpereira/projects/ntbuilder_cat/logs/req-%j.out
#SBATCH -e /mnt/home/mpereira/projects/ntbuilder_cat/logs/req-%j.err
set -u
source /etc/profile.d/modules.sh
module load python/3.12          # scipy/numpy/spglib; o python do sistema não tem
SITE=/mnt/home/mpereira/projects/ntbuilder_cat/site
export NTB_CATALOGUE_DB=${NTB_CATALOGUE_DB:-/mnt/home/mpereira/projects/ntbuilder_cat/catalogue.db}
cd "$SITE"
python3 request_job.py
