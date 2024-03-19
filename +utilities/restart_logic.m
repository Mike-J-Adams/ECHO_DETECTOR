function [iStart] = restart_logic(PATH2OUTPUT,PATH2DATA)
% Logic to restart script at last output 
%
% Last updated by Mike Adams
% 2024-03-19
% define output directory and check if there is existing output
  if isfolder(PATH2OUTPUT)
     % check for existing MAT files, get timestamp of the last one, and 
     % compare that with WAV file names to determine where to resume
     [~, matNames] = utilities.listFiles(PATH2OUTPUT, 'mat');
     [~,fileNames_wav] = utilities.listFiles(PATH2DATA,'wav');
     nFiles = numel(fileNames_wav);
     if ~isempty(matNames)
         lastTime = utilities.readDateTime(matNames{end});
         lastTime.Format = 'yyyyMMdd_HHmmss';
         timeStrParts = strsplit(char(lastTime),'_');
         leftOff = contains(fileNames_wav,timeStrParts{1}) & contains(fileNames_wav,timeStrParts{2});
         iStart = find(leftOff) + 1;
     else
         iStart = 1;
     end
     if iStart <= nFiles
         fprintf('Resuming at recording #%d/%d ("%s")\n',iStart,nFiles,fileNames_wav{iStart})
     else
         return
     end
  else
     % create new output folder and start from beginning
     mkdir(dirPath_out)
     iStart = 1;
     fprintf('Starting new extractions for deployment "%s"\n',depName) 
  end