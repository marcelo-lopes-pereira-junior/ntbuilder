"""
Envio de e-mail do site — pelo relay que o cluster já usa.

A máquina já manda e-mail: é assim que o Slurm avisa fim de job, via o Gmail
do laboratório (`lab.nanoeng@gmail.com`, só saída).  O site usa o mesmo
caminho, chamando ``sendmail``/``msmtp``, em vez de carregar uma senha de
aplicativo dentro do processo web.  Nenhuma credencial passa por aqui nem
pelo repositório.

Se o envio falhar — relay fora, permissão no ``msmtprc`` —, a mensagem é
gravada em ``web/data/mail/`` e o pedido segue em pé: o link continua válido e
a mensagem pode ser reenviada depois com

    cat web/data/mail/<arquivo>.eml | msmtp -t
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import time
from email.message import EmailMessage
from email.utils import formatdate
from pathlib import Path

_log = logging.getLogger("ntbuilder.mail")

MAIL_FROM = os.environ.get("NTB_MAIL_FROM", "lab.nanoeng@gmail.com")
MAIL_NAME = os.environ.get("NTB_MAIL_NAME", "NTBuilder — NanoEng")
SPOOL = Path(os.environ.get("NTB_MAIL_SPOOL",
                            Path(__file__).parent.parent / "data" / "mail"))
# Só desligado de propósito (desenvolvimento); vazio = ligado.
DISABLED = os.environ.get("NTB_MAIL_DISABLE", "").lower() in ("1", "yes", "true")


# /etc/msmtprc é root:slurm 0640 — legível pelo Slurm, que manda os avisos de
# job, e não pelo usuário do site.  Em vez de duplicar a senha de aplicativo,
# aponte esta variável para uma cópia do arquivo que o usuário do site possa
# ler (instalada pelo root, fora do repositório).
CONFIG = os.environ.get("NTB_MSMTP_CONFIG", "")


def _agent() -> list[str] | None:
    for cand in ("/usr/sbin/sendmail", "/usr/lib/sendmail"):
        if Path(cand).exists():
            return [cand, "-t", "-i"]
    msmtp = shutil.which("msmtp")
    if msmtp:
        return [msmtp, "-t"] + (["-C", CONFIG] if CONFIG else [])
    return None


def _spool(msg: EmailMessage, why: str) -> Path:
    SPOOL.mkdir(parents=True, exist_ok=True)
    name = f"{time.strftime('%Y%m%d-%H%M%S')}-{abs(hash(msg['Subject'])) % 10**6}.eml"
    path = SPOOL / name
    path.write_bytes(bytes(msg))
    _log.warning("e-mail não enviado (%s); guardado em %s", why, path)
    return path


def send(to: str, subject: str, body: str) -> tuple[bool, str]:
    """Manda a mensagem.  Devolve (enviou, motivo se não enviou)."""
    msg = EmailMessage()
    msg["From"] = f"{MAIL_NAME} <{MAIL_FROM}>"
    msg["To"] = to
    msg["Subject"] = subject
    msg["Date"] = formatdate(localtime=True)
    msg.set_content(body)

    if DISABLED:
        return False, "envio desligado (NTB_MAIL_DISABLE)"
    agent = _agent()
    if agent is None:
        _spool(msg, "nenhum sendmail/msmtp no sistema")
        return False, "nenhum agente de envio no sistema"
    try:
        r = subprocess.run(agent, input=bytes(msg), capture_output=True, timeout=30)
    except Exception as e:                       # timeout, agente travado
        _spool(msg, f"{type(e).__name__}: {e}")
        return False, f"{type(e).__name__}: {e}"
    if r.returncode != 0:
        why = (r.stderr or r.stdout or b"").decode(errors="replace").strip()
        _spool(msg, why or f"código {r.returncode}")
        return False, why or f"código {r.returncode}"
    return True, ""


def _num(n: int) -> str:
    """Milhar com ponto, como se escreve em português."""
    return f"{int(n):,}".replace(",", ".")


def _fmt_bytes(n: int) -> str:
    for unit in ("B", "kB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} GB"


def queued(to: str, name: str, protocol: str, est: dict, filter_text: str,
           base_url: str) -> tuple[bool, str]:
    who = f", {name}" if name else ""
    return send(to, f"[NTBuilder] pedido {protocol} recebido", f"""Olá{who},

seu pedido ao catálogo de nanotubos do NTBuilder foi recebido.

  Protocolo:  {protocol}
  Tubos:      {_num(est['n_tubes'])}
  Átomos:     {_num(est['n_atoms'])}
  Tamanho:    ~{_fmt_bytes(est['bytes_zip'])} compactado
  Estimativa: ~{est['seconds'] / 60:.0f} min de construção

Acompanhe em:
  {base_url}/catalogo?protocolo={protocol}

O filtro pedido:
{filter_text}

Quando terminar, você recebe outro e-mail com o link do arquivo.  O link fica
válido por 7 dias.

--
NTBuilder — Laboratório NanoEng, Universidade de Brasília
""")


def ready(to: str, name: str, protocol: str, n_files: int, size: int,
          expires: str, base_url: str) -> tuple[bool, str]:
    who = f", {name}" if name else ""
    return send(to, f"[NTBuilder] pedido {protocol} pronto", f"""Olá{who},

seu pedido está pronto.

  Protocolo: {protocol}
  Arquivos:  {_num(n_files)}
  Tamanho:   {_fmt_bytes(size)}
  Expira em: {expires[:10]}

Baixe em:
  {base_url}/api/cat/download/{protocol}

Dentro do arquivo vão também o manifest.csv, com a origem e os números de cada
tubo, e o CITATION.txt com as citações das bases de onde vieram as camadas —
por favor cite-as no trabalho.

--
NTBuilder — Laboratório NanoEng, Universidade de Brasília
""")


def failed(to: str, name: str, protocol: str, message: str,
           base_url: str) -> tuple[bool, str]:
    who = f", {name}" if name else ""
    return send(to, f"[NTBuilder] pedido {protocol} falhou", f"""Olá{who},

seu pedido {protocol} não pôde ser atendido:

  {message}

Se o filtro era grande, tente dividi-lo; se o erro parecer do nosso lado,
responda este e-mail com o protocolo e nós olhamos.

--
NTBuilder — Laboratório NanoEng, Universidade de Brasília
""")
