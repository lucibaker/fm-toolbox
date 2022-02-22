function [] = errorbar_CI(xdat,ydat,varargin)
%xdat,ydat of the form, [avg,CI];
%optional input, color in form of k,m,g,c, etc.

if ~isempty(varargin)
    col=varargin{1};

    if length(col)==3
        %col={'o-','color',col};
        col={'o','color',col};

     end

    if length(col)==1
        col=[col 'o'];
    end
else
    col='o';
end
    
if length(varargin)>1
    %opt=varargin(2:end);
    opt=varargin;
    col={};
else
    opt={};
end


[r,c]=size(xdat);
if c==1
    
    if isempty(col)
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),opt{:})%,'markerfacecolor',col)
    elseif iscell(col)
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),col{:},opt{:})%,'markerfacecolor',col)    
    else
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),col,opt{:})%,'markerfacecolor',col)
    end
else
    
    if isempty(col)
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),xdat(:,1)-xdat(:,2),xdat(:,3)-xdat(:,1),opt{:})%,'markersize',3)%,'markerfacecolor',col,'markersize',5)
    elseif iscell(col)
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),xdat(:,1)-xdat(:,2),xdat(:,3)-xdat(:,1),col{:},opt{:})%,'markersize',3)%,'markerfacecolor',col,'markersize',5)    
    else
        errorbar(xdat(:,1),ydat(:,1),ydat(:,1)-ydat(:,2),ydat(:,3)-ydat(:,1),xdat(:,1)-xdat(:,2),xdat(:,3)-xdat(:,1),col,opt{:})%,'markersize',3)%,'markerfacecolor',col,'markersize',5)
    end
end

end