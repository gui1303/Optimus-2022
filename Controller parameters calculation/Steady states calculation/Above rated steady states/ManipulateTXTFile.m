% Function: ManipulateTXTFile String replacement in .txt file
%
%
% -----------------------------
% Usage:
% -------------
% ManipulateTXTFile(TXTFile,StringToReplace,NewString)
% ------------
% Input:
% -------------
% TXTFile           String with name of .txt file
% StringToReplace   String which has to be replaced
% NewString         String with replacement
% Mode              (optional) String; If "regex", regexprep will be used
%                   for string replacement instead of strrep.
% ------------
% Output:
% ------------
% -
% ------------
% Needs:
% ------------
% -
% ------------
% Modified:
% -------------
% * Martin Koch on 21-01-2021
%   - added "regex" mode
% * David Schlipf (TTI GmbH) on 27-Apr-2017
%   - regexprep replaced by strrep
%   - replace Backslash removed, should be outside of this function
% * David Schlipf on 29-Dec-2011
%   - replace Backslash
% ------------
% ToDo:
% -------------
%
% -----------
% Created:
% David Schlipf on 05-Dec-2010
% (c) Universitaet Stuttgart and sowento (TTI-GmbH)
% ----------------------------------


function n = ManipulateTXTFile(TXTFile,StringToReplace,NewString,Mode)

useRegex = nargin == 4 && Mode == "regex";

TempTXTFile     = [TXTFile(1:end-4),'_temp',TXTFile(end-3:end)];

fid             = fopen(TXTFile);
fidTemp         = fopen(TempTXTFile,'w+');

n               = 0;

while ~feof(fid)
    s         = fgetl(fid);
    if useRegex
        sTemp = regexprep(s,StringToReplace,NewString);
    else
        sTemp = strrep(s,StringToReplace,NewString);
    end
    fprintf(fidTemp,'%s\r\n',sTemp);
    if ~strcmp(s,sTemp)
        n = n+1;
    end
    
end

fclose(fid);
fclose(fidTemp);
recycle         = 'on';
delete(TXTFile);
movefile(TempTXTFile,TXTFile);