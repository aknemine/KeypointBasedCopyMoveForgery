function deger = ortalama(dizi)

    [height, width] = size(ortalama);
    sayici = 0;
    toplam = 0;
    for i = 1 : height
        for j = 1 : width
            toplam = toplam + dizi(i,j);
            sayici = sayici + 1;
        end
    end
    deger = (toplam/sayici);
end