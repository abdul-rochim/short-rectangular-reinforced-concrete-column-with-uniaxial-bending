using Plots
using Random

#mutu material
fc = 40 #Mpa mutu beton
fy = 420 #Mpa mutu baja
Es = 200000 #Mpa Modulus baja tul
Ec = 4700*sqrt(fc)
eps_c = 0.003

#dimensi
c1 = 500 #mm (h)
c2 = c1 #mm (b)
cover = 40 #mm
y = c1/2 #pusat massa
Ag = c1 * c2

#baja tulangan
#asumsi single
dia = 19 #mm
n = 4*4 #buah, jumlah yg dipasang perimeter
sengkang = 13 #mm
As = n*0.25*pi*dia^2
d = c1 - cover - sengkang - 0.5*dia

#atur posisi tulangan
jumlah_baris = (n - 4)/4 + 2
y_top = cover + sengkang + 0.5*dia
#jarak_antar_tulangan = (c1 - 2*y_top)/(jumlah_baris-1)

y_posisi =  range(y_top, stop = c1-y_top, length = Int64(jumlah_baris))

#jumlah iterasi
n_iter = 1000
iter_at_Pn_pure = zeros(n_iter)
iter_at_Pn_pure_used = 0
iter_at_Pn = zeros(n_iter)
iter_at_Pn_used = 0

#data beban [Mu] [Pu]
Mu = rand(-300:10:300, 100) #kN.m
Pu = rand(-200:100:4000, 100) #kN

#analisa
c = range(0.0001*c1, c1*1.31, n_iter)
a = zeros(n_iter)
eps_s = zeros(Int64(jumlah_baris))
fs = zeros(Int64(jumlah_baris))
y_pos_to_NA = zeros(Int64(jumlah_baris))

As_per_baris = zeros(Int64(jumlah_baris))
As_negatif = zeros(Int64(jumlah_baris))
Pn_ = zeros(Int64(jumlah_baris))
Mn_ = zeros(Int64(jumlah_baris))

sum_As_negatif = 0

#Beton tertekan
Cc = zeros(Int64(n_iter))

Pn = zeros(Int64(n_iter))
Mn = zeros(Int64(n_iter))
phi_Pn = zeros(Int64(n_iter))
phi_Mn = zeros(Int64(n_iter))

#epsilon baja tulangan untuk mencari nilai phi
eps_t = zeros(Int64(n_iter))
max_eps_t = 0

#asumsi betha1
betha1 = 0.8
if fc >= 17 && fc <= 28
    betha1 = 0.85
elseif fc >= 55
    betha1 = 0.65
else
    betha1 = 0.85 - 0.05*(fc - 28)/7
end

#faktor reduksi kekuatan
phi = zeros(Int64(n_iter))

#regangan baja tul posisi serat tarik terluar
if fy == 420
    eps_ty = 0.002
else
    eps_ty = fy/Es
end

#iterasi
for i = 1 : n_iter
    a[i] = c[i] * betha1

    for j = 1 : Int64(jumlah_baris)
        y_pos_to_NA[j] = y_posisi[j] - c[i]
        eps_s[j] = eps_c * y_pos_to_NA[j]/c[i]
        eps_t[j] = eps_s[j]

        if y_posisi[j] > y_top && y_posisi[j] < c1-y_top
            As_per_baris[j] = As/n*2
        else
            As_per_baris[j] = As/n*((n-4)/4+2)
        end

        if y_posisi[j] < a[i]
            if y_posisi[j] > y_top && y_posisi[j] < c1-y_top
                As_negatif[j] = As/n*2
            else
                As_negatif[j] = As/n*((n-4)/4+2)
            end
        else
            #do nothing
        end

        if y_posisi[j] < c[i]
            if abs(eps_s[j]) >= fy/Es
                eps_s[j] = -fy/Es
            else
                #do nothing
            end
        else
            if eps_s[j] >= fy/Es
                eps_s[j] = fy/Es
            else
                #do nothing
            end
        end

        fs[j] = eps_s[j] * Es
        #if y_posisi[j] > c[i]
            Pn_[j] = - As_per_baris[j] * fs[j]
        #else
        #    Pn_[j] = - As_per_baris[j] * (fs[j] - 0.85* fc)
        #end

        Mn_[j] = As_per_baris[j] * fs[j] * (y_posisi[j] - y)
    end

    global sum_As_negatif = sum(As_negatif) #mm2

    if a[i] <= c1
        Cc[i] = 0.85*fc*(c2*a[i]-sum_As_negatif) #N
    else
        Cc[i] = 0.85*fc*(c2*c1-sum_As_negatif) #N
    end
    if a[i] >= c1
        if abs(Cc[i] - Cc[i-1]) < 0.0001
            Cc[i] = Cc[i-1]
        end
    end

    Pn[i] = (sum(Pn_) + Cc[i])/1000 #kN
    Mn[i] = (sum(Mn_) + Cc[i]*(y-a[i]/2))/1000000 #kN.m

#    println("maksimum regangan baja tul terluar : ", maximum(eps_t))
    global max_eps_t = abs(maximum(eps_t))

    if max_eps_t <= eps_ty
        phi[i] = 0.65
    elseif max_eps_t >= (eps_ty + 0.003)
        phi[i] = 0.9
    else
        phi[i] = 0.65 + 0.25*(max_eps_t - eps_ty)/0.003
    end

    phi_Pn[i] = phi[i] * Pn[i]
    phi_Mn[i] = phi[i] * Mn[i]
end

#pure compression
Pn_max = (0.85*fc*(Ag-sum_As_negatif) + fy*sum_As_negatif)/1000
Pn_pure = 0.65*Pn_max
println(Pn_max)
println(Pn_pure)
println(Pn_pure*0.8)
#phi_Pn_pure = phi[n_iter] * Pn_pure
#println("Pn_pure manual calc : ", Pn_pure)
#println("Pn_pure iterasi : ", Pn[n_iter])
#println("phi Pn pure manual calc : ", phi_Pn_pure)
#println()

#println("iter average at pure: ", Int64(average(iter_at_Pn_pure)))
#other pure compression
#iter_at_Pn_pure_used = Int64(maximum(iter_at_Pn_pure))
#println("iter maks at pure: ", iter_at_Pn_pure_used)
#Pn_pure_other = Pn[iter_at_Pn_pure_used] #(0.85*fc*(Ag*betha1-sum_As_negatif) + sum(Pn_))/1000
#phi_Pn_pure_other = phi[iter_at_Pn_pure_used] * Pn_pure_other
#println("Pn pure at iter : ", Pn_pure_other)
#println("phi Pn pure at iter : ", phi_Pn_pure_other)

#Plots
x_dummy = rand(-300:30:500, 100)
y_dummy = rand(-100:70:2000, 100)
p_dummy = plot(x_dummy, y_dummy, label = " ", color =:white)

p_pure_a = plot(p_dummy, [Mn[n_iter], 0], [Pn[n_iter], Pn_max], label = nothing, color =:red, style =:dot #=dash=#, lw = 2) #"Pn maks, Mn", color =:red, style =:dash, lw = 1)
p_pure_b = plot(p_pure_a, [-Mn[n_iter], 0], [Pn[n_iter], Pn_max], label = nothing, color =:red, style =:dot #=dash=#, lw = 2) #"Pn maks, -Mn", color =:red, style =:dash, lw = 1)

pb = scatter(p_pure_b, Mu, Pu, label = nothing, color =:blue, markersize = 1) #"Pu, Mu", color =:blue, markersize = 1)
p1 = plot(pb, Mn, Pn, label = nothing, color =:red, lw = 2, style =:dot) #"Pn, +Mn", color =:red, lw = 2, style =:dot)
p2 = plot(p1, -Mn, Pn, label = nothing, color =:red, lw = 2, style =:dot) #"Pn, -Mn",color =:red, lw = 2, style =:dot)

phi_p_pure_a = plot(p2, [phi_Mn[n_iter], 0], [phi_Pn[n_iter], Pn_pure], color =:orange, lw = 2, label = nothing) #"Pn pure, +phi Mn pure")
phi_p_pure_b = plot(phi_p_pure_a, [-phi_Mn[n_iter], 0], [phi_Pn[n_iter], Pn_pure], color =:orange, lw = 2, label = nothing) #"Pn pure, -phi Mn pure")

p3 = plot(phi_p_pure_b, phi_Mn, phi_Pn, label = nothing, color =:orange, lw = 2) #"phi_Pn, +phi_Mn", color =:orange, lw = 2)
plot!(p3, -phi_Mn, phi_Pn, title = "Diagram Interaksi", color =:orange, label = nothing, lw = 2) #"phi_Mn, -phi_Mn", lw = 2)
xlabel!("Momen (kN.m)")
ylabel!("Aksial (kN)")
