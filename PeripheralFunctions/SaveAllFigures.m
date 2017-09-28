function SaveAllFigures(basefn,filetype,maximize)

if nargin < 1
    basefn = 'Unknown'; % 'Unknown' filename
end
if nargin < 2
    filetype = 'eps'; % save as color eps by default
end
if nargin < 3
    maximize = 0; % do not maximize figures by default
end

ChildList = sort(get(0,'Children'));

for cnum = 1:length(ChildList)
    if strncmp(get(ChildList(cnum),'Type'),'figure',6)
        fig = ChildList(cnum);
        % maximize:
        if maximize
            units = get(fig,'units');
            set(fig,'units','normalized','outerposition',[0 0 1 1]);
            set(fig,'units',units);
        end
        
        % make sure outputs at same size as appears on screen -- crucial!!
        set(fig,'PaperPositionMode','auto');

        % set filename
        fn = [basefn '_FigureSave_' num2str(cnum) '.' filetype];
        
        % save and close
        if strcmp(filetype,'eps')
            print(fig,'-depsc2',fn)
%             saveas(fig, ['FigureSave_', fn, '_', num2str(fig), '.' filetype],'epsc'); % COLOR EPS
        else
            saveas(fig,fn);
        end
        close(fig)
    end
end

end