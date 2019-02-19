sumaS2 = Float64(0.36398547250893341852488170816398)
pi = Float64(pi)

# obliczam sumę częściową "od końca", aby uzyskać większą dokładność

#pierwsza metoda
function S2(k)
    sumaPozytywnych = 0.0
    sumaNegatywnych = 0.0
    for i = 1:k
        j = k+1-i
        if j % 2 == 1
            sumaPozytywnych += 1.0 / Float64(j^2+1)
        else
            sumaNegatywnych += 1.0 / Float64(j^2+1)
        end
    end
    return sumaPozytywnych - sumaNegatywnych
end

#druga metoda
function trocheUlepszonyS2(k)
    sumaEfektywna = Float64(pi)*Float64(pi) / 12.0
    sumaPozytywnych = 0.0
    sumaNegatywnych = 0.0
    for i = 1:k
        j = k+1-i
        if j % 2 == 1
            sumaPozytywnych += 1.0 / Float64(j*j*(j*j+1))
        else
            sumaNegatywnych += 1.0 / Float64(j*j*(j*j+1))
        end
    end
    return sumaEfektywna + sumaNegatywnych - sumaPozytywnych
end

#trzecia metoda
function bardziejUlepszonyS2(k)
    sumaEfektywna = pi*pi / 12.0 - pi*pi*pi*pi*7 / 720.0
    sumaPozytywnych = 0.0
    sumaNegatywnych = 0.0
    for i = 1:k
        j = k+1-i
        if j % 2 == 1
            sumaPozytywnych += 1.0 / Float64(j*j*j*j*(j*j+1))
        else
            sumaNegatywnych += 1.0 / Float64(j*j*j*j*(j*j+1))
        end
    end
    return sumaEfektywna + sumaPozytywnych - sumaNegatywnych
end

function statystykiSumCzesciowychS2prec10(A, fun)
    for k in A
        suma = fun(k)
        @printf("Iter: %d \t%.10f\t%e\n",k,suma,abs(sumaS2-suma))
    end
end

function statystykiSumCzesciowychS2prec16(A, fun)
    for k in A
        suma = fun(k)
        @printf("Iter: %d \t%.16f\t%e\n",k,suma,abs(sumaS2-suma))
    end
end

# funkcja szukająca k takiego, że a_{k+1} <= e, gdzie e to błąd dokładności
# to znaczy, jaką sumę częściową trzeba policzyć, by uzyskać zadaną dokładność
function szukajSumyCzesciowej(j=0, n=2, prec=10)
    k = 1
    bound = Float64(10.0^Float64(-prec)/2.0)
    while 1.0 / Float64(k^j*(k^n+1)) > bound
        k+=1
    end
    return k-1
end

function statystykiBledowBezwzglS2(A,j,n)
    for k in A
        @printf("Liczba kroków: %d, \t błąd bezwzgl.: %e\n", k, 1.0 / Float64((k+1)^j*((k+1)^n+1)))
    end
end
