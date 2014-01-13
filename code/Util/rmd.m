function v = rmd(w)
v(1) = w(1);
c = 2;
for i=2:length(w)
    if w(i) > w(i-1)
        v(c) = w(i);
        c = c + 1;
    end
end