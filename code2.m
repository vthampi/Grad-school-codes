%cdata = [2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 5.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 4.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 4.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 4.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;2.000 2.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 3.000 1.000;1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000]; %K-W
%W-K

%Make your matrix ^ Write it in reverse (top is first entry
%so make sure to account for this in code. Easier to load MAT)!!

%cdata2 = [2 3 3 3 3 3 3 3 3 3 1;2 3 3 3 3 3 3 3 3 3 1;2 3 3 3 3 3 3 3 3 3 1;2 3 3 3 3 3 3 3 3 3 1;2 4 3 3 3 3 3 3 3 3 1;2 4 4 4 7 7 7 7 7 7 1;2 4 4 4 7 7 7 7 7 7 1;2 2 2 6 6 6 6 6 6 6 1;2 2 2 6 6 6 6 6 6 6 1;2 2 2 6 6 6 6 6 6 6 6;2 2 2 6 6 6 6 6 6 6 6];
%W-RHO

%cdata3 = [1 3 3 3 3 3 3 3 3 3 3;1 3 3 3 3 3 3 3 3 3 3;1 3 3 3 3 3 3 3 3 3 3;1 3 3 3 3 3 3 3 3 3 3;1 5 5 5 5 5 5 5 5 5 5;1 4 4 4 4 4 4 4 4 4 4;1 4 4 4 4 4 4 4 4 4 4;8 2 2 2 2 2 6 6 6 6 6;2 2 2 2 2 2 6 6 6 6 6;2 2 2 2 2 2 6 6 6 6 6;2 2 2 2 2 2 6 6 6 6 6];
%K-RHO

cdata4 = [1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1;2 3 3 3 3 3 3 3 3 3 1;2 4 3 3 3 3 3 3 3 3 1;2 4 4 7 7 7 7 7 7 7 1;2 4 4 7 7 7 7 7 7 7 1;2 2 6 6 6 6 6 6 6 6 1;2 9 6 6 6 6 6 6 6 6 6];
%I-W

x = 0:0.1:1;
%What does x go from?
y = 5:-0.5:0;
%What Does y go from? Remember the reverse thing above? Do it again.
clf
h = heatmap(x, y, cdata4,'Title','I-W Parameter Plane','XLabel','W','YLabel','I','Colormap',jet(10),'CellLabelColor','none');
%Colormap tells you what color style you're using. There are many, check
%online (literally, Google "Colormap Matlab").
%CellLabelColor is the thing that writes down the value in the cell itself.
%Type none to turn it off.
% everything is pretty self explanatory from here.
