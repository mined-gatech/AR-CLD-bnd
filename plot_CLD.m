function plot_CLD(xbin,f,varargin)

    flip = 0;
    clim = 0;
    set_theme = 0;
    if nargin > 2
        for ii=1:length(varargin)
            if strcmp(varargin{ii},'flip') == 1
                flip = 1;
            end
            if strcmp(varargin{ii},'clim') == 1
                clim = 1;
                lims = varargin{ii+1};
            end
            if strcmp(varargin{ii},'theme') == 1
                set_theme = 1;
                theme = varargin{ii+1};
            end
        end
    end

    ang_step = 2*pi/(size(f,1)-1);
        
    [r,t] = meshgrid(xbin,0:ang_step:2*pi);
    x = r.*cos(t);
    y = r.*sin(t);
    
    [~,h] = contourf(x,y,f,128); axis image; axis off
    set(h,'edgecolor','none');
    if set_theme && exist('brewermap','file')
        map = brewermap(128,theme);
        if flip 
            map = flipud(map);
        elseif strcmp(theme,'Spectral') || strcmp(theme,'RdBu')
            fprintf('Warning: Spectral and RdBu usually need to be flipped\n');
        end
        colormap(map)
    end
    
    % set limits
    if clim
        caxis(lims);
    end
end
    