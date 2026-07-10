#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

# ==============================================================================
# INSTALLER STANDALONE PER deseq_auto_import.R
#
# Installa tutto ciò che serve allo script R:
#   - R
#   - DESeq2
#   - stringr
#   - pheatmap
#   - ggplot2
#   - matrixStats
#   - RColorBrewer
#   - dplyr
#
# Il pacchetto "tools" è incluso nell'installazione base di R e non deve essere
# installato separatamente.
#
# Ambiente creato:
#   ~/miniforge3/envs/nanopore_deseq2
#
# Esecuzione consigliata:
#   chmod +x install_deseq2_R_env.sh
#   ./install_deseq2_R_env.sh
#
# Lo script accetta anche:
#   sudo ./install_deseq2_R_env.sh
#
# Se viene usato sudo, l'ambiente viene comunque installato nella home
# dell'utente che ha invocato sudo, non nella home di root.
#
# IMPORTANTE:
#   se nanopore_deseq2 esiste già, viene eliminato e ricreato da zero.
# ==============================================================================

ENV_NAME="nanopore_deseq2"

R_VERSION="4.5"
DESEQ2_VERSION="1.50.2"

MINIFORGE_VERSION="26.3.2-3"

log() {
    printf '\n\033[1;34m[%s]\033[0m %s\n' "$(date '+%F %T')" "$*"
}

ok() {
    printf '\033[1;32m✅ %s\033[0m\n' "$*"
}

warn() {
    printf '\033[1;33m⚠️  %s\033[0m\n' "$*" >&2
}

die() {
    printf '\033[1;31m❌ ERRORE: %s\033[0m\n' "$*" >&2
    exit 1
}

on_error() {
    local status=$?
    local line="${1:-?}"

    printf '\n\033[1;31m❌ Installazione interrotta alla riga %s, codice %s.\033[0m\n' \
        "${line}" "${status}" >&2

    exit "${status}"
}

trap 'on_error "${LINENO}"' ERR

# ==============================================================================
# 1. UTENTE E ARCHITETTURA
# ==============================================================================

if [[ "${EUID}" -eq 0 ]]; then
    TARGET_USER="${SUDO_USER:-}"

    if [[ -z "${TARGET_USER}" || "${TARGET_USER}" == "root" ]]; then
        die "Avvia lo script come utente normale oppure tramite sudo da un utente normale."
    fi
else
    TARGET_USER="${USER}"
fi

command -v getent >/dev/null 2>&1 \
    || die "Il comando getent non è disponibile."

TARGET_HOME="$(getent passwd "${TARGET_USER}" | cut -d: -f6)"
TARGET_GROUP="$(id -gn "${TARGET_USER}")"

[[ -n "${TARGET_HOME}" && -d "${TARGET_HOME}" ]] \
    || die "Home non valida per ${TARGET_USER}: ${TARGET_HOME}"

case "$(uname -m)" in
    x86_64)
        ARCH="x86_64"
        MINIFORGE_SHA256="848194851a98903134187fbb4ab50efe87b003e0c0f808f97644b7524a62bf2c"
        ;;
    aarch64)
        ARCH="aarch64"
        MINIFORGE_SHA256="2c113a69297e612b01ca0f320c22a3107a11f2ab9b573d79ac868a175945ce29"
        ;;
    *)
        die "Architettura non supportata: $(uname -m)"
        ;;
esac

run_user() {
    if [[ "${EUID}" -eq 0 ]]; then
        sudo -H -u "${TARGET_USER}" env \
            HOME="${TARGET_HOME}" \
            USER="${TARGET_USER}" \
            LOGNAME="${TARGET_USER}" \
            "$@"
    else
        env \
            HOME="${TARGET_HOME}" \
            USER="${TARGET_USER}" \
            LOGNAME="${TARGET_USER}" \
            "$@"
    fi
}

run_apt() {
    if [[ "${EUID}" -eq 0 ]]; then
        apt-get "$@"
    else
        sudo apt-get "$@"
    fi
}

fix_owner() {
    if [[ "${EUID}" -eq 0 ]]; then
        chown -R "${TARGET_USER}:${TARGET_GROUP}" "$@"
    fi
}

# ==============================================================================
# 2. PERCORSI
# ==============================================================================

MINIFORGE_DIR="${TARGET_HOME}/miniforge3"
CONDA="${MINIFORGE_DIR}/bin/conda"
MAMBA="${MINIFORGE_DIR}/bin/mamba"

ENV_DIR="${MINIFORGE_DIR}/envs/${ENV_NAME}"

INSTALL_DIR="${TARGET_HOME}/nanopore_R_install"
LOCK_DIR="${INSTALL_DIR}/locks"
YAML_FILE="${INSTALL_DIR}/${ENV_NAME}.yml"

MINIFORGE_INSTALLER="Miniforge3-${MINIFORGE_VERSION}-Linux-${ARCH}.sh"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_INSTALLER}"

log "Utente destinatario : ${TARGET_USER}"
log "Home utente         : ${TARGET_HOME}"
log "Architettura        : ${ARCH}"
log "Miniforge           : ${MINIFORGE_DIR}"
log "Ambiente R          : ${ENV_DIR}"

# ==============================================================================
# 3. DIPENDENZE DEL SISTEMA
# ==============================================================================

command -v apt-get >/dev/null 2>&1 \
    || die "Questo installer richiede Ubuntu o Debian con apt-get."

log "Installazione delle dipendenze di sistema..."

export DEBIAN_FRONTEND=noninteractive

run_apt update

run_apt install -y --no-install-recommends \
    bash \
    bzip2 \
    ca-certificates \
    coreutils \
    curl \
    fontconfig \
    fonts-dejavu-core \
    gzip \
    tar \
    xz-utils

ok "Dipendenze di sistema disponibili."

# ==============================================================================
# 4. MINIFORGE
# ==============================================================================

if [[ -x "${CONDA}" && -x "${MAMBA}" ]]; then
    log "Miniforge è già installato e verrà riutilizzato."
else
    if [[ -e "${MINIFORGE_DIR}" ]]; then
        BACKUP_DIR="${MINIFORGE_DIR}.incompleto.$(date '+%Y%m%d_%H%M%S')"

        warn "${MINIFORGE_DIR} esiste ma non è un'installazione completa."
        warn "La directory verrà spostata in ${BACKUP_DIR}."

        mv "${MINIFORGE_DIR}" "${BACKUP_DIR}"
        fix_owner "${BACKUP_DIR}"
    fi

    TMP_DIR="$(mktemp -d)"
    trap 'rm -rf "${TMP_DIR}"' EXIT

    INSTALLER_PATH="${TMP_DIR}/${MINIFORGE_INSTALLER}"

    log "Download di Miniforge ${MINIFORGE_VERSION}..."

    curl \
        --fail \
        --location \
        --retry 5 \
        --retry-delay 3 \
        --output "${INSTALLER_PATH}" \
        "${MINIFORGE_URL}"

    [[ -s "${INSTALLER_PATH}" ]] \
        || die "L'installer Miniforge scaricato è vuoto."

    log "Verifica dello SHA-256 ufficiale..."

    printf '%s  %s\n' \
        "${MINIFORGE_SHA256}" \
        "${INSTALLER_PATH}" \
        | sha256sum --check -

    log "Installazione di Miniforge..."

    run_user bash "${INSTALLER_PATH}" \
        -b \
        -p "${MINIFORGE_DIR}"

    rm -rf "${TMP_DIR}"
    trap - EXIT
fi

[[ -x "${CONDA}" ]] \
    || die "Conda non trovato: ${CONDA}"

[[ -x "${MAMBA}" ]] \
    || die "Mamba non trovato: ${MAMBA}"

fix_owner "${MINIFORGE_DIR}"

log "Versioni del gestore di ambienti:"
run_user "${CONDA}" --version
run_user "${MAMBA}" --version

ok "Miniforge pronto."

# ==============================================================================
# 5. DEFINIZIONE RIPRODUCIBILE DELL'AMBIENTE
#
# I pacchetti vengono risolti tutti insieme in un ambiente pulito.
# Non viene eseguito alcun aggiornamento incrementale.
# ==============================================================================

mkdir -p "${INSTALL_DIR}" "${LOCK_DIR}"
fix_owner "${INSTALL_DIR}"

cat > "${YAML_FILE}" <<EOF
name: ${ENV_NAME}

channels:
  - conda-forge
  - bioconda
  - nodefaults

dependencies:
  - r-base=${R_VERSION}
  - bioconductor-deseq2=${DESEQ2_VERSION}
  - r-stringr
  - r-pheatmap
  - r-ggplot2
  - r-matrixstats
  - r-rcolorbrewer
  - r-dplyr
EOF

fix_owner "${YAML_FILE}"

ok "Definizione YAML creata: ${YAML_FILE}"

# ==============================================================================
# 6. RIMOZIONE E RICREAZIONE DELL'AMBIENTE
# ==============================================================================

if [[ -e "${ENV_DIR}" ]]; then
    log "Rimozione dell'ambiente precedente ${ENV_NAME}..."

    if ! run_user "${MAMBA}" env remove \
        --prefix "${ENV_DIR}" \
        --yes; then

        warn "La rimozione Conda non è riuscita; elimino la directory residua."
        rm -rf "${ENV_DIR}"
    fi
fi

log "Creazione pulita dell'ambiente ${ENV_NAME}..."

run_user "${MAMBA}" env create \
    --prefix "${ENV_DIR}" \
    --file "${YAML_FILE}" \
    --yes

[[ -x "${ENV_DIR}/bin/Rscript" ]] \
    || die "Rscript non trovato dopo la creazione dell'ambiente."

fix_owner "${ENV_DIR}"

ok "Ambiente ${ENV_NAME} creato."

# ==============================================================================
# 7. VERIFICA DEI PACCHETTI USATI DA deseq_auto_import.R
# ==============================================================================

log "Verifica dei pacchetti R..."

run_user env \
    PATH="${ENV_DIR}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    R_HOME="${ENV_DIR}/lib/R" \
    "${ENV_DIR}/bin/Rscript" --vanilla -e '
required <- c(
    "DESeq2",
    "stringr",
    "tools",
    "pheatmap",
    "ggplot2",
    "matrixStats",
    "RColorBrewer",
    "dplyr"
)

missing <- required[
    !vapply(
        required,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing) > 0L) {
    stop(
        "Pacchetti R mancanti: ",
        paste(missing, collapse = ", ")
    )
}

cat("R: ", R.version.string, "\n", sep = "")

for (package in required) {
    cat(
        package,
        ": ",
        as.character(packageVersion(package)),
        "\n",
        sep = ""
    )
}
'

ok "Tutti i pacchetti richiesti dallo script R sono disponibili."

# ==============================================================================
# 8. FILE DI LOCK
# ==============================================================================

log "Esportazione dei file di lock..."

run_user "${CONDA}" list \
    --prefix "${ENV_DIR}" \
    --explicit \
    > "${LOCK_DIR}/${ENV_NAME}-${ARCH}.lock.txt"

run_user "${CONDA}" env export \
    --prefix "${ENV_DIR}" \
    > "${LOCK_DIR}/${ENV_NAME}-full.yml"

fix_owner "${INSTALL_DIR}"

ok "File di lock salvati in ${LOCK_DIR}."

# ==============================================================================
# 9. ISTRUZIONI FINALI
# ==============================================================================

cat <<EOF

==============================================================================
✅ INSTALLAZIONE R/DESEQ2 COMPLETATA
==============================================================================

Ambiente:
    ${ENV_DIR}

Attivazione manuale:

    source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
    conda activate "${ENV_DIR}"

Esecuzione manuale dello script R:

    Rscript \\
        /home/${TARGET_USER}/reference/deseq_auto_import.R \\
        /percorso/cartella/07_deseq2 \\
        PB \\
        SANO

La pipeline Bash attiva automaticamente questo ambiente prima di eseguire R.

Definizione YAML:
    ${YAML_FILE}

File di lock:
    ${LOCK_DIR}

==============================================================================
EOF
