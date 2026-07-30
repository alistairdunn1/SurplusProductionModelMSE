@echo off
rem =====================================
rem Build and Check SurplusProductionModelMSE R Package (Windows)
rem =====================================

rem Check if R is in PATH
where R >nul 2>nul
if errorlevel 1 (
  echo R is not in your PATH. Please install R and add it to your PATH.
  exit /b 1
)

rem Generate documentation with roxygen2
echo Generating documentation...
rm -rf SurplusProductionModelMSE.Rcheck
call R --vanilla < run-roxygen.R
if errorlevel 1 exit /b 1

rem Build the package (includes vignettes)
echo Building package with vignettes...
call R CMD build --force .
if errorlevel 1 exit /b 1

rem Install the package from the built tar.gz
echo Installing package...
for %%f in (SurplusProductionModelMSE_*.tar.gz) do (
  call R CMD INSTALL "%%f"
  if errorlevel 1 exit /b 1
)

rem Check the built package
echo Checking package...
for %%f in (SurplusProductionModelMSE_*.tar.gz) do (
  call R CMD check "%%f"
  if errorlevel 1 exit /b 1
)

echo.
echo =====================================
echo Build completed successfully!
echo =====================================
