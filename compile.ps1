$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location (Join-Path $Root "Paper")

if (-not (Test-Path ".\ghandali_index2.tex")) {
  Write-Error "ghandali_index2.tex not found in Paper\. Flatten Paper\PAPER if needed."
}

Get-ChildItem -ErrorAction SilentlyContinue *.aux,*.bbl,*.blg,*.log,*.out,*.toc,*.synctex.gz |
  Remove-Item -Force -ErrorAction SilentlyContinue

pdflatex -interaction=nonstopmode ghandali_index2.tex
if ($LASTEXITCODE -ne 0) { Write-Warning "pdflatex pass 1 had errors (see main.log)" }

bibtex ghandali_index2
pdflatex -interaction=nonstopmode ghandali_index2.tex
pdflatex -interaction=nonstopmode ghandali_index2.tex

if (Test-Path ".\ghandali_index2.pdf") {
  $pdf = Get-Item ".\ghandali_index2.pdf"
  Write-Host "OK: $($pdf.FullName)  ($($pdf.Length) bytes)"
} else {
  Write-Error "PDF not produced. Check Paper\main.log"
}
