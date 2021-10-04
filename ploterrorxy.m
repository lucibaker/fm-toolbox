function[] = ploterrorxy(x,y,errx,erry)
%this function plots a scatter plot with x and y error bars. 

errorbar(x,y,erry,'.k-','MarkerSize',10)
hold on
herrorbar(x,y,errx,'.k-')
end
