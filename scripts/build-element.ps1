$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$compilerDirectory = 'C:\msys64\mingw64\bin'
if (-not (Test-Path -LiteralPath "$compilerDirectory\g++.exe")) {
    throw "C++ compiler not found in $compilerDirectory"
}
$env:PATH = "$compilerDirectory;$env:PATH"
& cmake -S $projectRoot -B "$projectRoot\build" -G 'MinGW Makefiles' "-DCMAKE_CXX_COMPILER=$compilerDirectory\g++.exe" "-DCMAKE_MAKE_PROGRAM=$compilerDirectory\mingw32-make.exe" -DCMAKE_BUILD_TYPE=Release
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed.' }
& cmake --build "$projectRoot\build" --target element_demo test_element test_mesh_reader --parallel 2
if ($LASTEXITCODE -ne 0) { throw 'Compilation failed.' }
& ctest --test-dir "$projectRoot\build" -R '^(element|mesh_reader)$' --output-on-failure
if ($LASTEXITCODE -ne 0) { throw 'Element/reader tests failed.' }
Write-Output 'Element and reader checks passed. See docs/element-hand-check.md before viewing reference answers.'
