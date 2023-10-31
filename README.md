# Millig
C++ legacy project about math modeling of milling process 
Требования:
- Поддержка OpenMP
- наличие установленной boost (libboost-all-dev)

Компиляция: 
```
g++ -std=c++14 -fopenmp diplom.cpp -o diplom
```
Запуск:
./diplom 

Опциональные аргументы запуска:

* -L [dbl] - феноменологический коэффициент
* -typePAV [int] - тип пов. активного вещ-ва  Al2O3drinding: -1 , Dry grinding: 0, SAS: 1, Isopropyl: 2 , Ethanol: 3
* -Con [dbl] масса SiC г.
* -oborot [dbl] 250 или 400 оборотов в минуту
* -densBalls [dbl] плотность шаров кг/м^3
* -densParticle [dbl] плотность частиц кг/м^3
* -typeMill [int] тип мельницы, 1 если модель Иванникова, 0 - модель Бабкина
* -P [dbl] прединтегральный коэффициент
* -massBallsAndParticles [dbl] [dbl] масса шаров и масса частиц
* -massRatio [dbl] отношение массы шаров и массы частиц
* -sizeBall [dbl] размер мелющих шаров, мм
* -searchSize [dbl] размер который необходимо получить (в микронах)
* -avStart [dbl] average particle size (23,7)micron
* -Tmax [dbl] время дробления в секундах
