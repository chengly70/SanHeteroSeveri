function b=bm(v)        % Beta for activation gating variable m in Na current
b=8000*exp(-0.056*(v+66));
end