"""
Lightweight SQLite layer for NTBuilder Web.

Stores:
  - downloads : job_id, email, name, timestamp
"""
from __future__ import annotations

import asyncio
import sqlite3
from datetime import datetime
from pathlib import Path

_DB_PATH = Path(__file__).parent.parent / "data" / "ntbuilder.db"


def _get_conn() -> sqlite3.Connection:
    _DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(_DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


async def init_db() -> None:
    """Create tables on first startup."""
    loop = asyncio.get_event_loop()
    await loop.run_in_executor(None, _init_db_sync)


def _init_db_sync() -> None:
    with _get_conn() as conn:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS downloads (
                id        INTEGER PRIMARY KEY AUTOINCREMENT,
                job_id    TEXT NOT NULL,
                email     TEXT NOT NULL,
                name      TEXT,
                timestamp TEXT NOT NULL
            )
        """)
        conn.commit()


async def save_registration(job_id: str, email: str, name: str = "") -> None:
    loop = asyncio.get_event_loop()
    await loop.run_in_executor(None, _save_sync, job_id, email, name)


def _save_sync(job_id: str, email: str, name: str) -> None:
    ts = datetime.utcnow().isoformat()
    with _get_conn() as conn:
        conn.execute(
            "INSERT INTO downloads (job_id, email, name, timestamp) VALUES (?,?,?,?)",
            (job_id, email, name, ts),
        )
        conn.commit()


async def get_job_email(job_id: str) -> str | None:
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(None, _get_email_sync, job_id)


def _get_email_sync(job_id: str) -> str | None:
    with _get_conn() as conn:
        row = conn.execute(
            "SELECT email FROM downloads WHERE job_id=? LIMIT 1", (job_id,)
        ).fetchone()
    return row["email"] if row else None


async def get_all_registrations() -> list[dict]:
    """Admin helper: return all download registrations."""
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(None, _get_all_sync)


def _get_all_sync() -> list[dict]:
    with _get_conn() as conn:
        rows = conn.execute(
            "SELECT job_id, email, name, timestamp FROM downloads ORDER BY timestamp DESC"
        ).fetchall()
    return [dict(r) for r in rows]
