function weights = initialize_Glorot_V2(numIn, numOut,className)

arguments
    numIn
    numOut
    className = 'single'
end

Z = 2*rand([numIn,numOut],className) - 1;
bound = sqrt(6 / (numIn + numOut));

weights = bound * Z;
weights = dlarray(weights);

end