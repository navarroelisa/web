param([string]$Pandoc = $env:PANDOC)
$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
  if (-not $Pandoc) {
    $command = Get-Command pandoc -ErrorAction SilentlyContinue
    if ($command) { $Pandoc = $command.Source }
    else {
      $bundled = Join-Path $env:ProgramFiles 'RStudio\resources\app\bin\quarto\bin\tools\pandoc.exe'
      if (Test-Path $bundled) { $Pandoc = $bundled }
      else { throw 'Pandoc not found. Use ./build.ps1 -Pandoc <path-to-pandoc> or render with rmarkdown::render_site().' }
    }
  }
  New-Item -ItemType Directory -Force docs | Out-Null
  foreach ($page in @('index','elisa','noemi','labmembers','papers','news','outreach','research')) {
    & $Pandoc "$page.Rmd" --from=markdown --to=html5 --standalone --template=site-template.html --output="docs/$page.html"
    if ($LASTEXITCODE -ne 0) { throw "Pandoc failed for $page" }
  }
  foreach ($folder in @('assets','images')) {
    New-Item -ItemType Directory -Force "docs/$folder" | Out-Null
    Copy-Item "$folder/*" "docs/$folder" -Recurse -Force
  }
  Write-Host 'Built eight pages in docs/.'
} finally { Pop-Location }
