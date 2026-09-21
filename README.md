# MitoGlia Lab website

Scientific website for Elisa Navarro and Noemí Esteras, Universidad Complutense de Madrid.

## Editing and building

Edit the `.Rmd` files. The shared layout is `site-template.html`, the responsive warm BuPu theme is `assets/site.css`, and `assets/site.js` handles the accessible mobile navigation. The font uses the visitor's system sans-serif (Segoe UI on Windows); no font service is required.

From PowerShell, run `./build.ps1`. It uses Pandoc from PATH or the installed RStudio bundle. An explicit binary can be supplied with `./build.ps1 -Pandoc <path>`. The pages contain Markdown and HTML, with no executable R chunks, so this build does not require R packages. R Markdown users can also use `rmarkdown::render_site()` with the shared template configured in `_site.yml`.

The generated site is in `docs/`. Preview it with a local HTTP server rooted at `docs/`. Keep source changes and generated HTML/assets together when publishing using the repository's existing workflow.

## Publications

The 21 September 2026 update uses the public ORCID records for Elisa (0000-0002-4056-7146) and Noemí (0000-0002-7938-6131), supplemented by Crossref DOI metadata. A concise source snapshot is in `data/publications-2026-09-21.json`. The publication list is static, so visitors do not depend on ORCID being available. Update `papers.Rmd` and rebuild when new works appear. DOI duplicates are collapsed; verified journal versions replace earlier preprints. Remaining preprints are explicitly labelled.

## Team

Current PhD students: Mamen, Lucía, María and David. Santiago is listed under “PhDs defended in the lab”; Sara Carmona is in Alumni. No defence date or current affiliation was inferred.
