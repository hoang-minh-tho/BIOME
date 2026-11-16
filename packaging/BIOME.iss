#define AppName "BIOME"
#define AppVersion "1.0.0"
#define AppPublisher "BIOME"
#define AppExeName "BIOME.exe"
#define AppSource AddBackslash(SourcePath) + "..\dist\BIOME"

[Setup]
AppId={{7E6BFAC8-4F68-44A9-9A2A-58154458F1BD}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
DefaultDirName={localappdata}\Programs\BIOME
DefaultGroupName={#AppName}
PrivilegesRequired=lowest
ArchitecturesInstallIn64BitMode=x64
WizardStyle=modern
SetupIconFile=..\assets\biome.ico
OutputDir=..\dist
OutputBaseFilename=BIOME-Setup
Compression=lzma2
SolidCompression=yes
UninstallDisplayIcon={app}\{#AppExeName}

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"

[InstallDelete]
Type: files; Name: "{app}\Databio16S.exe"
Type: files; Name: "{autodesktop}\Databio Amplicon.lnk"
Type: files; Name: "{userprograms}\Databio Amplicon\Databio Amplicon.lnk"

[Files]
Source: "{#AppSource}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\BIOME"; Filename: "{app}\{#AppExeName}"
Name: "{autodesktop}\BIOME"; Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; Description: "Launch BIOME"; Flags: postinstall nowait skipifsilent
