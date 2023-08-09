function  clean_up

%
%   clean up for tidy experiment ending...
%
%
%   hot tip: put a shortcut to this in the MATLAB toolbar for easy access
%   when something goes wrong (and you have two displays)...
%
%
%   ella wufong - March 2014 - microfish@fishmonkey.com.au
%

Priority(0); % reset priority

RestrictKeysForKbCheck([]); % allow KbCheck to listen to all keys again...

ListenChar(1);  %   re-enable keyboard input

Screen('CloseAll');     % close Psychtoolbox window

clear PsychImaging;     % clear PsychImaging, e.g. for VideoSwitcher

ShowCursor;  % restore mouse cursor

%PsychPortAudio('Close');    % close all open audio devices

commandwindow;  % return focus to MATLAB command window


end

