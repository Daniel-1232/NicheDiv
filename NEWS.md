# NicheDiv 1.1.0 (2026-10-21)

## Accessible-area delineation

- Added user-defined accessible-area polygons to crop.background.buffered().
- Added buffer.method = "custom" and the custom.polygon argument.
- Added support for sf, sfc, terra SpatVector, and spatial polygon files.
- Added automatic polygon CRS transformation and geometry validation.
- Supported multipolygons, disconnected regions, and polygon holes.
- Custom polygons are used directly without additional buffering.
- Preserved the existing hull, points, alpha, and bbox methods.
- Added permanent regression tests for custom polygon handling.

# NicheDiv 1.0.0 (2026-10-21)

## First stable release

- First stable release of the NicheDiv R package for multivariate niche divergence analyses using discriminant analysis of principal components (DAPC).
- Updated ClimateNAr integration to version 3.2.0.
- Added cross-platform ClimateNAr installation and extraction support for Windows, Linux, and macOS.
- Added GitHub Actions validation for Windows, Ubuntu Linux, macOS Intel, and macOS Apple Silicon.
- Improved package portability and cross-platform testing.

