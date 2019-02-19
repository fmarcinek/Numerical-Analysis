# WĘZŁY:

# Funkcja obliczajaca zera n-tego wielomianu Czebyszewa przeskalowanego do przedziału [a,b]
function wezlyCzebyszewa(n,a,b)
    wezly = Array{BigFloat}(n)

    for k in n:-1:1
         wezly[n+1-k] = 0.5*(a+b) + 0.5*(b-a)*cos(((2.0*k-1.0)/(2.0*n))*BigFloat(pi))
    end

    return wezly
end

# Funkcja obliczajaca n wezlow rownoodleglych dla przedziału [a,b] i n >= 2
function wezlyRownoodlegle(n,a,b)
    wezly = Array{BigFloat}(n)
    h = (b-a) / (n-1)

    for k in 1:n
         wezly[k] = a + (k-1) * h
    end

    return wezly
end

# Funkcja zwracajaca wartosci funkcji f w podanych wezlach (args)
function podajWartosci(f, args)
    return map(f, args)
end


# INTERPOLACJA WIELOMIANOWA:

# Funkcja wyliczajaca ilorazy różnicowe dla funkcji zadanej
# tablica argumentow (args) oraz tablica wartosci w punktach args (values)
function ilorazRoznicowy(args, values)
    n = length(args)

    diffValues = copy(values)

    for j in 1:(n-1)
        for i in n:(-1):(1+j)
            diffValues[i] = (diffValues[i] - diffValues[i-1]) / (args[i] - args[i-j])
        end
    end

    return diffValues;
end


# Funkcja interpolujaca funkcje w wezlach (args), ktora to funkcja zadana jest tablica
# argumentow (args) i odpowiadajacych im wartosci w tablicy (values)
function interpolacjaNewtona(args, values)
    b = ilorazRoznicowy(args, values)   # wspolczynniki wielomianu interpolujacego w postaci Newtona
    n = length(args)

    w = Array{Function}(n)      # tablica wielomianow (n-ty wielomian bedzie wielomianem interpolujacym)
    w[1] = _ -> b[1]

    p = Array{Function}(n)
    p[1] = _ -> 1.0

    for i in 2:n
        p[i] = x -> p[i-1](x) * (x - args[i-1])
        w[i] = x -> w[i-1](x) + b[i] * p[i](x)
    end

    return w[n]            # wielomian interpolujacy
end


# INTERPOLACJA FUNKCJAMI SKLEJANYMI:

# Funkcja obliczajaca drugie pochodne w wezlach naturalnej funkcji sklejanej
# trzeciego stopnia dla wezlow danych w tablicy args i wartosci interpolowanej
# funkcji w tych wezlach w tablicy values.
# Algorytm poznany na wykladzie, posiadajacy zlozonosc liniowa.
function drugiePochodne(args, values)
    n = length(args)
    u = Array{BigFloat}(n-1)
    q = Array{BigFloat}(n-1)

    u[1] = 0.0
    q[1] = 0.0

    for k in 1:(n-2)
        λ_k = (args[k+1]-args[k])/(args[k+2]-args[k+1]) + (args[k+1]-args[k])
        d_k = 6.0 * ilorazRoznicowy(args, [args[k],args[k+1],args[k+2]])[3]
        p = λ_k * q[k] + 2.0
        q[k+1] = (λ_k - 1.0) / p
        u[k+1] = (d_k - λ_k * u[k]) / p
    end

    M = Array{BigFloat}(n)   # tablica obliczanych pochodnych
    M[1] = 0.0
    M[n] = 0.0
    M[n-1] = u[n-1]

    for k in (n-3):(-1):1
        M[k+1] = u[k+1] + q[k+1]*M[k+2]
    end

    return M
end


# Funkcja obliczajaca wartosc w punkcie c naturalnej funkcji sklejanej trzeciego
# stopnia dla funkcji f, ktora jest interpolowana w wezlach zapisanych w tablicy args.
# args - tablica wezlow
# values - tablica wartosci funkcji f
# M - tablica drugich pochodnych wspomnianej funkcji sklejanej w wezlach args
function auxiliarySpline(c, args, values, M)

    k = 1
    while values[k+1] < c
        k += 1
    end

    h = values[k+1] - values[k]

    return ( (M[k]*(values[k+1] - c)^3.0) / 6.0 + (M[k+1]*(c - values[k])^3.0) / 6.0 +
                (args[k]-(M[k]*(h^2.0)) / 6.0)*(values[k+1] - c) +
                (args[k+1]-(M[k+1]*(h^2.0)) / 6.0)*(c - values[k]) ) / h
end


# INTERPOLACJA FUNKCJI ODWROTNEJ:

# Mozliwe tryby:
CZ = "Czebyszew"
RO = "Rownoodl"
# Funkcja interpolujaca metoda Newtona wartosc funkcji odwrotnej do f na przedziale [a,b] w punkcie c
# w n wezlach Czebyszewa (tryb="Czebyszew") lub w n wezlach rownoodleglych (tryb="Rownoodl")
# dla n >= 2
function interpolacjaOdwrotna(f, c, a, b, n, tryb)
    args = Array{BigFloat}(n)
    values = Array{BigFloat}(n)

    if tryb == CZ
        values = wezlyCzebyszewa(n, a, b)
    elseif tryb == RO
        values = wezlyRownoodlegle(n, a, b)
    else
        return
    end

    args = podajWartosci(f, values)

    w = interpolacjaNewtona(args, values)   # wielomian interpolacyjny
    return w(c)
end

# Funkcja obliczajaca blad pomiedzy interpolowana wartoscia a rzeczywista
function funkcjaBłędu(f, c, a, b, n, tryb, wart)
    wart = BigFloat(wart)
    return abs(BigFloat(interpolacjaOdwrotna(f, c, a, b, n, tryb)) - wart)
end

# Funkcja rysujaca wykres bledu miedzy wartoscia interpolowana a rzeczywista dla sposobu doboru wezlow w zaleznosci od ich liczby
function errorPlot(f, c, a, b, n0, n1, tryb, wart)
    plot_args = n0:n1
    n = length(plot_args)

    function auxErrorFunc(x)
        return funkcjaBłędu(f, c, a, b, x, tryb, wart)
    end

    func_values = Array{BigFloat}(n)
    func_values = map(auxErrorFunc, plot_args)

    func = scatter(;x=plot_args, y=func_values, mode="lines", name="funkcja interpolowana")

    if tryb == CZ
        layout = Layout(;title="Interpolacja w węzłach Czebyszewa",xaxis=attr(title="liczba węzłów Czebyszewa"),yaxis=attr(title="wartości błędu bezwzgl."))
    else
        layout = Layout(;title="Interpolacja w węzłach równoodległych",xaxis=attr(title="liczba węzłów równoodległych"),yaxis=attr(title="wartości błędu bezwzgl."))
    end

    plot(func,layout)
end
