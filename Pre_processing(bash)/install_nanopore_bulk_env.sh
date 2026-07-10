#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

# ==============================================================================
# INSTALLER RIPRODUCIBILE — NANOPORE BULK RNA-SEQ
#
# Esecuzione supportata:
#
#   ./install_nanopore_bulk_env.sh
#
# oppure:
#
#   sudo ./install_nanopore_bulk_env.sh
#
# Se viene usato sudo, Miniforge e gli ambienti vengono comunque installati
# nella home dell'utente che ha invocato sudo, non nella home di root.
#
# Lo script:
#   1. installa le dipendenze Ubuntu/Debian;
#   2. installa Miniforge 26.3.2-3 con SHA-256 ufficiale incorporato;
#   3. elimina e ricrea da zero i tre ambienti della pipeline;
#   4. verifica i programmi con i comandi ufficiali;
#   5. esporta file di lock riproducibili.
#
# Ambienti:
#   nanopore_bulk_stable
#       NanoPlot, Pychopper, minimap2, samtools, featureCounts
#
#   nanopore_fastcat
#       fastcat isolato, secondo la documentazione Oxford Nanopore
#
#   nanopore_deseq2
#       R, DESeq2 e pacchetti R
#
# La pipeline deve attivare solamente nanopore_bulk_stable.
# ==============================================================================

MAIN_ENV_NAME="nanopore_bulk_stable"
FASTCAT_ENV_NAME="nanopore_fastcat"
R_ENV_NAME="nanopore_deseq2"

PYTHON_VERSION="3.11"
PYSAM_VERSION="0.24.0"
NANOPLOT_VERSION="1.47.1"
PYCHOPPER_VERSION="2.7.10"
MINIMAP2_VERSION="2.31"
SAMTOOLS_VERSION="1.23.1"
SUBREAD_VERSION="2.1.1"

FASTCAT_VERSION="1.0.1"

R_VERSION="4.5"
DESEQ2_VERSION="1.50.2"
APEGLM_VERSION="1.32.0"
ENHANCEDVOLCANO_VERSION="1.28.2"

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
# 1. UTENTE DESTINATARIO E ARCHITETTURA
# ==============================================================================

if [[ "${EUID}" -eq 0 ]]; then
    TARGET_USER="${SUDO_USER:-}"

    if [[ -z "${TARGET_USER}" || "${TARGET_USER}" == "root" ]]; then
        die "Esegui lo script come utente normale oppure mediante sudo da un utente normale."
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

MAIN_ENV_DIR="${MINIFORGE_DIR}/envs/${MAIN_ENV_NAME}"
FASTCAT_ENV_DIR="${MINIFORGE_DIR}/envs/${FASTCAT_ENV_NAME}"
R_ENV_DIR="${MINIFORGE_DIR}/envs/${R_ENV_NAME}"

INSTALL_DIR="${TARGET_HOME}/nanopore_bulk_install"
LOCK_DIR="${INSTALL_DIR}/locks"

MINIFORGE_INSTALLER="Miniforge3-${MINIFORGE_VERSION}-Linux-${ARCH}.sh"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_INSTALLER}"

log "Utente destinatario : ${TARGET_USER}"
log "Home utente         : ${TARGET_HOME}"
log "Architettura        : ${ARCH}"
log "Miniforge           : ${MINIFORGE_DIR}"
log "Ambiente principale : ${MAIN_ENV_NAME}"
log "Ambiente fastcat    : ${FASTCAT_ENV_NAME}"
log "Ambiente R/DESeq2   : ${R_ENV_NAME}"

# ==============================================================================
# 3. DIPENDENZE UBUNTU/DEBIAN
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
    gzip \
    tar \
    xz-utils

ok "Dipendenze di sistema disponibili."

# ==============================================================================
# 4. INSTALLAZIONE DI MINIFORGE
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

[[ -x "${CONDA}" ]] || die "Conda non trovato in ${CONDA}."
[[ -x "${MAMBA}" ]] || die "Mamba non trovato in ${MAMBA}."

fix_owner "${MINIFORGE_DIR}"

log "Versioni del gestore di ambienti:"
run_user "${CONDA}" --version
run_user "${MAMBA}" --version

ok "Miniforge pronto."

# ==============================================================================
# 5. RIMOZIONE PULITA DEI VECCHI AMBIENTI
# ==============================================================================

remove_environment() {
    local env_name="$1"
    local env_dir="$2"

    if [[ ! -e "${env_dir}" ]]; then
        return 0
    fi

    log "Rimozione dell'ambiente precedente ${env_name}..."

    if ! run_user "${MAMBA}" env remove \
        --prefix "${env_dir}" \
        --yes; then

        warn "Rimozione Conda non riuscita; elimino la directory residua."
        rm -rf "${env_dir}"
    fi
}

remove_environment "${MAIN_ENV_NAME}" "${MAIN_ENV_DIR}"
remove_environment "${FASTCAT_ENV_NAME}" "${FASTCAT_ENV_DIR}"
remove_environment "${R_ENV_NAME}" "${R_ENV_DIR}"

rm -rf \
    "${TARGET_HOME}/nanopore_bulk_conda" \
    "${TARGET_HOME}/.local/lib/nanopore_bulk"

mkdir -p "${LOCK_DIR}"
fix_owner "${INSTALL_DIR}"

# ==============================================================================
# 6. AMBIENTE PRINCIPALE
# ==============================================================================

log "Creazione dell'ambiente principale..."

run_user "${MAMBA}" create \
    --yes \
    --prefix "${MAIN_ENV_DIR}" \
    --override-channels \
    --strict-channel-priority \
    --channel conda-forge \
    --channel bioconda \
    "python=${PYTHON_VERSION}" \
    "pysam=${PYSAM_VERSION}" \
    "nanoplot=${NANOPLOT_VERSION}" \
    "pychopper=${PYCHOPPER_VERSION}" \
    "minimap2=${MINIMAP2_VERSION}" \
    "samtools=${SAMTOOLS_VERSION}" \
    "subread=${SUBREAD_VERSION}" \
    gzip \
    pigz

[[ -x "${MAIN_ENV_DIR}/bin/python" ]] \
    || die "L'ambiente principale non è stato creato correttamente."

ok "Ambiente principale creato."

# ==============================================================================
# 7. AMBIENTE FASTCAT
# ==============================================================================

log "Creazione dell'ambiente fastcat..."

run_user "${MAMBA}" create \
    --yes \
    --prefix "${FASTCAT_ENV_DIR}" \
    --override-channels \
    --strict-channel-priority \
    --channel conda-forge \
    --channel bioconda \
    --channel nanoporetech \
    "nanoporetech::fastcat=${FASTCAT_VERSION}"

[[ -x "${FASTCAT_ENV_DIR}/bin/fastcat" ]] \
    || die "L'ambiente fastcat non è stato creato correttamente."

ok "Ambiente fastcat creato."

# ==============================================================================
# 8. AMBIENTE R / DESEQ2
# ==============================================================================

log "Creazione dell'ambiente R/DESeq2..."

run_user "${MAMBA}" create \
    --yes \
    --prefix "${R_ENV_DIR}" \
    --override-channels \
    --strict-channel-priority \
    --channel conda-forge \
    --channel bioconda \
    "r-base=${R_VERSION}" \
    "bioconductor-deseq2=${DESEQ2_VERSION}" \
    bioconductor-biocparallel \
    "bioconductor-apeglm=${APEGLM_VERSION}" \
    "bioconductor-enhancedvolcano=${ENHANCEDVOLCANO_VERSION}" \
    r-ashr \
    r-data.table \
    r-dplyr \
    r-tidyr \
    r-tibble \
    r-readr \
    r-stringr \
    r-purrr \
    r-ggplot2 \
    r-ggrepel \
    r-pheatmap \
    r-rcolorbrewer \
    r-openxlsx \
    r-optparse \
    r-jsonlite \
    r-cowplot \
    r-patchwork \
    r-scales

[[ -x "${R_ENV_DIR}/bin/Rscript" ]] \
    || die "L'ambiente R/DESeq2 non è stato creato correttamente."

ok "Ambiente R/DESeq2 creato."

# ==============================================================================
# 9. WRAPPER NELL'AMBIENTE PRINCIPALE
# ==============================================================================

log "Collegamento di fastcat, R e Rscript all'ambiente principale..."

cat > "${MAIN_ENV_DIR}/bin/fastcat" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${FASTCAT_ENV_DIR}/bin/fastcat" "\$@"
EOF

cat > "${MAIN_ENV_DIR}/bin/R" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="${R_ENV_DIR}/bin:\${PATH}"
export R_HOME="${R_ENV_DIR}/lib/R"
exec "${R_ENV_DIR}/bin/R" "\$@"
EOF

cat > "${MAIN_ENV_DIR}/bin/Rscript" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="${R_ENV_DIR}/bin:\${PATH}"
export R_HOME="${R_ENV_DIR}/lib/R"
exec "${R_ENV_DIR}/bin/Rscript" "\$@"
EOF

chmod 0755 \
    "${MAIN_ENV_DIR}/bin/fastcat" \
    "${MAIN_ENV_DIR}/bin/R" \
    "${MAIN_ENV_DIR}/bin/Rscript"

fix_owner \
    "${MAIN_ENV_DIR}/bin/fastcat" \
    "${MAIN_ENV_DIR}/bin/R" \
    "${MAIN_ENV_DIR}/bin/Rscript"

ok "Wrapper creati."

# ==============================================================================
# 10. CONTROLLI STANDARD
# ==============================================================================

log "Controllo dei programmi installati..."

run_user "${MAIN_ENV_DIR}/bin/python" --version
run_user "${MAIN_ENV_DIR}/bin/NanoPlot" --version
run_user "${MAIN_ENV_DIR}/bin/pychopper" --help >/dev/null
run_user "${MAIN_ENV_DIR}/bin/minimap2" --version
run_user "${MAIN_ENV_DIR}/bin/samtools" --version
run_user "${MAIN_ENV_DIR}/bin/featureCounts" -v
run_user "${MAIN_ENV_DIR}/bin/fastcat" fastq --version
run_user "${MAIN_ENV_DIR}/bin/Rscript" --version

run_user "${MAIN_ENV_DIR}/bin/Rscript" -e '
required <- c(
    "DESeq2",
    "BiocParallel",
    "apeglm",
    "ashr",
    "EnhancedVolcano",
    "data.table",
    "dplyr",
    "tidyr",
    "tibble",
    "readr",
    "stringr",
    "purrr",
    "ggplot2",
    "ggrepel",
    "pheatmap",
    "RColorBrewer",
    "openxlsx",
    "optparse",
    "jsonlite",
    "cowplot",
    "patchwork",
    "scales"
)

missing <- required[
    !vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing) > 0L) {
    stop(
        "Pacchetti R mancanti: ",
        paste(missing, collapse = ", ")
    )
}

cat(
    "DESeq2 ",
    as.character(packageVersion("DESeq2")),
    "\n",
    sep = ""
)
'

ok "Tutti i programmi richiesti sono disponibili."

# ==============================================================================
# 11. FILE DI LOCK
# ==============================================================================

log "Esportazione delle build installate..."

run_user "${CONDA}" list \
    --prefix "${MAIN_ENV_DIR}" \
    --explicit \
    > "${LOCK_DIR}/${MAIN_ENV_NAME}-${ARCH}.lock.txt"

run_user "${CONDA}" list \
    --prefix "${FASTCAT_ENV_DIR}" \
    --explicit \
    > "${LOCK_DIR}/${FASTCAT_ENV_NAME}-${ARCH}.lock.txt"

run_user "${CONDA}" list \
    --prefix "${R_ENV_DIR}" \
    --explicit \
    > "${LOCK_DIR}/${R_ENV_NAME}-${ARCH}.lock.txt"

fix_owner "${INSTALL_DIR}"

ok "File di lock creati in ${LOCK_DIR}."

# ==============================================================================
# 12. RISULTATO
# ==============================================================================

cat <<EOF

==============================================================================
INSTALLAZIONE COMPLETATA
==============================================================================

Miniforge:
    ${MINIFORGE_DIR}

Ambiente da attivare nella pipeline:
    ${MAIN_ENV_NAME}

Configurazione della pipeline:

    CONDA_BASE="${MINIFORGE_DIR}"
    source "\${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${MAIN_ENV_NAME}"

Sintassi fastcat:

    fastcat fastq \\
        --min-length 500 \\
        --min-qscore 10 \\
        input.fastq \\
        > output.fastq

File di lock:
    ${LOCK_DIR}

==============================================================================
EOF
