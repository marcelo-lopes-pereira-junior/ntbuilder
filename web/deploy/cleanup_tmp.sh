#!/bin/bash
# Cron job: delete job directories older than 24 hours.
# Schedule with: crontab -e
#   0 * * * *  /opt/ntbuilder/web/deploy/cleanup_tmp.sh >> /var/log/ntbuilder-cleanup.log 2>&1

TMP_DIR="/opt/ntbuilder/web/tmp"
MAX_AGE_HOURS=24

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Cleaning $TMP_DIR (older than ${MAX_AGE_HOURS}h)"

# Remove job subdirectories older than MAX_AGE_HOURS
find "$TMP_DIR" -mindepth 1 -maxdepth 1 -type d \
     -mmin "+$((MAX_AGE_HOURS * 60))" \
     -exec rm -rf {} \; -print

# Remove orphan .cif uploads older than 2 hours
find "$TMP_DIR" -maxdepth 1 -name "*.cif" \
     -mmin "+120" \
     -exec rm -f {} \; -print

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Done."
