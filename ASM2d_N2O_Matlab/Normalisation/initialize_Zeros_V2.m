function parameter = initialize_Zeros_V2(numState,className)

arguments
    numState
    className = 'single'
end

parameter = zeros([1, numState],className);
parameter = dlarray(parameter);

end