# compile.ps1 — build Paper/ghandali_index2.pdf from the repository root
# Usage:  .\compile.ps1
# Requires: pdflatex, bibtex on PATH (MiKTeX or TeX Live)

$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Paper = Join-Path $Root "Paper"
if (-not (Test-Path -LiteralPath $Paper)) {
    Write-Error "Paper\ directory not found next to compile.ps1 (expected: $Paper)"
}
Set-Location -LiteralPath $Paper

if (-not (Test-Path -LiteralPath ".\ghandali_index2.tex")) {
    Write-Error "ghandali_index2.tex not found in Paper\."
}
if (-not (Test-Path -LiteralPath ".\ghandali_index2.bib")) {
    Write-Error "ghandali_index2.bib not found in Paper\ (required for BibTeX)."
}

# Clean previous auxiliaries so this is a reproducible full rebuild
Get-ChildItem -ErrorAction SilentlyContinue `
    *.aux, *.bbl, *.blg, *.log, *.out, *.toc, *.synctex.gz, *.fls, *.fdb_latexmk |
    Remove-Item -Force -ErrorAction SilentlyContinue

function Invoke-TeX {
    param(
        [Parameter(Mandatory = $true)][string]$Exe,
        [Parameter(Mandatory = $true)][string[]]$TeXArgs
    )
    # Do not name a parameter $Args — it collides with PowerShell's automatic $Args.
    & $Exe @TeXArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Error "$Exe failed with exit code $LASTEXITCODE (see Paper\ghandali_index2.log / .blg)"
    }
}

# Order: pdflatex -> bibtex -> pdflatex -> pdflatex
Invoke-TeX -Exe "pdflatex" -TeXArgs @("-interaction=nonstopmode", "ghandali_index2.tex")
Invoke-TeX -Exe "bibtex"   -TeXArgs @("ghandali_index2")
Invoke-TeX -Exe "pdflatex" -TeXArgs @("-interaction=nonstopmode", "ghandali_index2.tex")
Invoke-TeX -Exe "pdflatex" -TeXArgs @("-interaction=nonstopmode", "ghandali_index2.tex")

$pdf = Join-Path $Paper "ghandali_index2.pdf"
if (-not (Test-Path -LiteralPath $pdf)) {
    Write-Error "PDF not produced. Inspect Paper\ghandali_index2.log"
}

$item = Get-Item -LiteralPath $pdf
Write-Host ("OK: {0}  ({1} bytes)" -f $item.FullName, $item.Length)
