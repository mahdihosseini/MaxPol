function img(i)

i = (double(i));
imagesc(i)
axis image
colormap gray
set(gca, 'Xtick', [], 'Ytick', [])