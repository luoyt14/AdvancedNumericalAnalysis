function H = Hilbert(dimension)
% create Hilbert Matrix
% input: the demension of the Hilbert Matrix
H = zeros(dimension, dimension);
for k=1:dimension
    for m=1:dimension
        H(k,m)=1/(k+m-1);
    end
end

end
