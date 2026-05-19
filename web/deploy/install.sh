#!/bin/bash
# ─── NTBuilder Web — Server Installation Script ───────────────────────────────
# Run as root on the server: bash install.sh
# Assumes Ubuntu/Debian with Python 3.10+, nginx already installed.
set -euo pipefail

REPO_DIR="/opt/ntbuilder"
VENV_DIR="$REPO_DIR/venv"
WEB_DIR="$REPO_DIR/web"
SERVICE_FILE="$WEB_DIR/deploy/ntbuilder-web.service"
NGINX_SNIPPET="$WEB_DIR/deploy/nginx-snippet.conf"

echo "=== NTBuilder Web Installer ==="

# 1. Create venv and install deps
echo "[1/5] Creating Python virtual environment..."
python3 -m venv "$VENV_DIR"
"$VENV_DIR/bin/pip" install --upgrade pip wheel
"$VENV_DIR/bin/pip" install -r "$WEB_DIR/requirements.txt"

# 2. Create tmp and data dirs with correct permissions
echo "[2/5] Creating data directories..."
mkdir -p "$WEB_DIR/tmp" "$WEB_DIR/data"
chown -R www-data:www-data "$WEB_DIR/tmp" "$WEB_DIR/data"
chmod 750 "$WEB_DIR/tmp" "$WEB_DIR/data"

# 3. Install systemd service
echo "[3/5] Installing systemd service..."
# Update paths in service file if needed
sed "s|/opt/ntbuilder|$REPO_DIR|g" "$SERVICE_FILE" \
    > /etc/systemd/system/ntbuilder-web.service
systemctl daemon-reload
systemctl enable ntbuilder-web
systemctl restart ntbuilder-web
systemctl status ntbuilder-web --no-pager

# 4. Install cron for tmp cleanup
echo "[4/5] Installing cleanup cron job..."
chmod +x "$WEB_DIR/deploy/cleanup_tmp.sh"
CRON_LINE="0 * * * * $WEB_DIR/deploy/cleanup_tmp.sh >> /var/log/ntbuilder-cleanup.log 2>&1"
(crontab -l 2>/dev/null | grep -v "cleanup_tmp"; echo "$CRON_LINE") | crontab -

# 5. Nginx hint
echo "[5/5] nginx configuration:"
echo "  Add the contents of $NGINX_SNIPPET to your server block, then:"
echo "  sudo nginx -t && sudo systemctl reload nginx"
echo ""
echo "=== Done! NTBuilder Web is running at http://127.0.0.1:8765 ==="
echo "    After configuring nginx: https://www.nanoeng.unb.br/ntbuilder"
