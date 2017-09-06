% Algoritmo recursivo para la serie de Fibonacci 
function y = fibonacci(x)
if x == 0
    y = 0;
elseif x == 1
    y = 1;
else
    y = fibonacci(x-1)+fibonacci(x-2);
end
end