::=========================================================================
:: By : Charlotte.HonG
@Echo Off
Title CHG_MakeFile V3.2.1
::=========================================================================
::�M�צW��
set project=Sift
::�D�{�����|(�i�d��)
set folder=%OneDrive%\Git Repository\%project%
::=========================================================================
::�ɮצW�� - �h�ɽsĶ(�۰ʸɤW".cpp")
set main=%project%_main
set file0=
set file1=SingleFile.cpp
set file2=
set file3=
set file4=
set file5=
set file6=
set file7=
set file9=
::=========================================================================
::�i���ﶵ
::C++�з�(�i�d�Źw�])
set std=c++17
::ĵ�i(�i�d��)
set Wall=-Wall
::�u�ƿﶵ[0-3](�i�d��or���[��L���O)
set O=-O2
::=========================================================================
::���Y�ɸ��|
set Inc0=%USERPROFILE%\Desktop\%project%\head
set Inc1=
set Inc2=
set Inc3=
set Inc4=
::�禡�w���|
set Ldir0=
set Ldir1=
set Ldir2=
set Ldir3=
set Ldir4=
::�禡�w
set lib0=
set lib1=
set lib2=
set lib3=
set lib4=
set lib5=
set lib6=
set lib7=
set lib9=
::=========================================================================
::�ץ��Ѽ�
if "%main%" EQU "" set main="%main%"
if "%folder%" EQU "" set folder=%cd%
if "%std%" NEQ "" set std=-std=%std%
if "%file0%" NEQ "" set file0="%file0%"
if "%file1%" NEQ "" set file1="%file1%"
if "%file2%" NEQ "" set file2="%file2%"
if "%file3%" NEQ "" set file3="%file3%"
if "%file4%" NEQ "" set file4="%file4%"
if "%file5%" NEQ "" set file5="%file5%"
if "%file6%" NEQ "" set file6="%file6%"
if "%file7%" NEQ "" set file7="%file7%"
if "%file8%" NEQ "" set file8="%file8%"
if "%file9%" NEQ "" set file9="%file9%"
if "%Inc0%" NEQ "" set Inc0=-I"%Inc0%"
if "%Inc1%" NEQ "" set Inc1=-I"%Inc1%"
if "%Inc2%" NEQ "" set Inc2=-I"%Inc2%"
if "%Inc3%" NEQ "" set Inc3=-I"%Inc3%"
if "%Inc4%" NEQ "" set Inc4=-I"%Inc4%"
if "%Ldir0%" NEQ "" set Ldir0=-L"%Ldir0%"
if "%Ldir1%" NEQ "" set Ldir1=-L"%Ldir1%"
if "%Ldir2%" NEQ "" set Ldir2=-L"%Ldir2%"
if "%Ldir3%" NEQ "" set Ldir3=-L"%Ldir3%"
if "%Ldir4%" NEQ "" set Ldir4=-L"%Ldir4%"
if "%lib0%" NEQ "" set lib0=-l"%lib0%"
if "%lib1%" NEQ "" set lib1=-l"%lib1%"
if "%lib2%" NEQ "" set lib2=-l"%lib2%"
if "%lib3%" NEQ "" set lib3=-l"%lib3%"
if "%lib4%" NEQ "" set lib4=-l"%lib4%"
if "%lib5%" NEQ "" set lib5=-l"%lib5%"
if "%lib6%" NEQ "" set lib6=-l"%lib6%"
if "%lib7%" NEQ "" set lib7=-l"%lib7%"
if "%lib8%" NEQ "" set lib8=-l"%lib8%"
if "%lib9%" NEQ "" set lib9=-l"%lib9%"
::=========================================================================
::�P�w�O�_���u�@�ؿ�
if "%folder%" == "%cd%" (
::�sĶ�ð���{��
g++ %Wall% %std% %O% %main%.cpp^
 %Inc0% %Inc1% %Inc2% %Inc3% %Inc4%^
 %Ldir0% %Ldir1% %Ldir2% %Ldir3% %Ldir4%^
 %file0% %file1% %file2% %file3% %file4%^
 %file5% %file6% %file7% %file8% %file9%^
 %lib0% %lib1% %lib2% %lib3% %lib4%^
 %lib5% %lib6% %lib7% %lib8% %lib9%^
 -o %main% && cmd /c "%main% & "
) else (
::�����u�@�ؿ�
cd "%folder%" && cmd /c "%~n0%~x0"
)
::=========================================================================
