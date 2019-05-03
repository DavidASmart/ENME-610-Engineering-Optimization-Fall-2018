function make_animated_gif(varargin)
%% Description:
%   function to build an animated gif
%% Instructions:
%   There a few different inputs needed 
%   at various points within a function to successfully create a gif. 
%
%   First, use make_animated_gif('clear') to start a new gif. 
%
%   Then for each frame of the gif, use make_animated_gif('snap', [handle]), 
%
%   and when done with all the frames, 
%   use make_animated_gif('write', filename, [delayTime], [loopCount]) 
%   to generate the gif file.
%% Inputs:
%   filename        - name of the file you want to make
%   delayTime       - Delay time between images 		 (s, default is 0)
%   loopCount       - Number of times to loop gif 		(default is infinite)
%% Outputs:
% 	[filename].gif	- the animated gif file

%%


global ANIMATED_GIF_IM ANIMATED_GIF_MAP ANIMATED_GIF_PREV

command = varargin{1};



if strcmp(command,'clear')
% ======================================================================= %
% ======================================================================= %
    ANIMATED_GIF_IM=[];
    ANIMATED_GIF_MAP=[];
% ======================================================================= %
% ======================================================================= %

    
    
elseif strcmp(command,'snap')
% ======================================================================= %
% ======================================================================= %
    if nargin>=2
        handle=varargin{2};
    else
        handle=gcf;
    end
    
    set(gcf,'Renderer','zbuffer') 
    
    if isempty(ANIMATED_GIF_IM)
        drawnow
        f = getframe(handle);
        [ANIMATED_GIF_IM,ANIMATED_GIF_MAP] = ...
        rgb2ind(f.cdata,256,'nodither');
    else
        f = getframe(handle);
        try
            [new_im,new_map] = rgb2ind(f.cdata,256,'nodither');             % Account for any new colors to be added to the map
            ANIMATED_GIF_MAP = [ANIMATED_GIF_MAP; 
            setdiff(new_map,ANIMATED_GIF_MAP,'rows')];
            ANIMATED_GIF_IM(:,:,1,end+1) = ...
            rgb2ind(f.cdata,ANIMATED_GIF_MAP,'nodither');
        catch
            disp('Skipping differently sized image')
        end
    end
% ======================================================================= %
% ======================================================================= %
    
    

elseif strcmp(command,'write')
% ======================================================================= %
% ======================================================================= %
    filename=varargin{2};
    if nargin>=3
        delayTime=varargin{3};
    end
    if nargin>=4
        loopCount=varargin{4};
    end
    
    [pathstr name ext]=fileparts(filename);
    if ~isdir(pathstr) && ~isempty(pathstr)
        error('Path does not exist: %s\n',pathstr);
    end

    disp(['Writing to: ' fullfile(pathstr,[name '.gif'])])
    
    if isempty(ANIMATED_GIF_IM) && isfield(ANIMATED_GIF_PREV,'IM')
        ANIMATED_GIF_IM = ANIMATED_GIF_PREV.IM;
        ANIMATED_GIF_MAP= ANIMATED_GIF_PREV.MAP;
    end
    im=ANIMATED_GIF_IM;
    map=ANIMATED_GIF_MAP;
    
    if ~exist('delayTime')
        delayTime=0;
    end
    if ~exist('loopCount')
        loopCount=inf;
    end
    
    if isempty(im)
        error('Nothing to write.')
    else
        imwrite(im,map,fullfile(pathstr,[name '.gif']),...
        'DelayTime',delayTime,'LoopCount',loopCount);
    end     
    ANIMATED_GIF_PREV.IM = ANIMATED_GIF_IM;
    ANIMATED_GIF_PREV.MAP= ANIMATED_GIF_MAP;
    ANIMATED_GIF_IM=[];
    ANIMATED_GIF_MAP=[];
% ======================================================================= %
% ======================================================================= %
    
    
else
% ======================================================================= %
% ======================================================================= %
    error('Command not understood: %s',command)
% ======================================================================= %
% ======================================================================= %
end

end
    