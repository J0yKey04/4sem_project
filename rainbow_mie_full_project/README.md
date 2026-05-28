# Rainbow Mie Full Project

Проект считает **реальное Mie scattering** для сферической капли воды и строит полную азимутальную радугу на 360° вокруг антисолнечной точки.

## Что именно моделируется

- капля воды считается идеальной сферой;
- показатель преломления воды зависит от длины волны;
- для каждой длины волны считаются коэффициенты Ми `a_n`, `b_n`;
- считаются амплитудные функции рассеяния `S1(theta)`, `S2(theta)`;
- интенсивность берётся как

```text
I(theta) = 1/2 ( |S1(theta)|^2 + |S2(theta)|^2 )
```

- угол картинки задаётся через

```text
alpha = 180° - theta
```

где `alpha` — угловое расстояние от антисолнечной точки. Основная радуга находится примерно около `alpha ≈ 42°`.

## Структура

```text
rainbow_mie_full_project/
├── CMakeLists.txt
├── include/
│   ├── dispersion.hpp
│   ├── image.hpp
│   └── mie.hpp
├── src/
│   ├── dispersion.cpp
│   ├── image.cpp
│   ├── main.cpp
│   └── mie.cpp
├── scripts/
│   └── plot_results.py
└── results/
```

## Запуск на macOS / Linux

Из папки проекта:

```bash
cmake -S . -B build
cmake --build build
./build/rainbow_mie
```

После запуска появятся:

```text
results/angular_intensity.csv
results/full_360_rainbow.ppm
```

Открыть PPM на macOS можно через Preview или конвертировать:

```bash
sips -s format png results/full_360_rainbow.ppm --out results/full_360_rainbow.png
```

## Построить графики

```bash
python3 scripts/plot_results.py
```

Будут созданы:

```text
results/rainbow_region_plot.png
results/full_backscatter_plot.png
```

Если matplotlib не установлен:

```bash
python3 -m pip install matplotlib
```

## Параметры запуска

```bash
./build/rainbow_mie [radius_um] [angle_samples] [width] [height]
```

Примеры:

```bash
./build/rainbow_mie 10 3601 1200 1200
./build/rainbow_mie 50 3601 1200 1200
./build/rainbow_mie 100 7201 1600 1600
```

## Как интерпретировать картинку

- центр картинки — антисолнечная точка;
- радиус от центра — угол `alpha`;
- полный круг — это 360° по азимуту;
- светлое кольцо около 40–43° — основная радуга;
- более слабая структура около 50–54° — область вторичной радуги;
- тонкие осцилляции — волновая интерференционная структура Mie scattering, то есть то, чего нет в чистой геометрической оптике.

## Важные ограничения

Это не трассировка отдельных лучей в OpenGL-сцене, а более фундаментальная угловая модель рассеяния одной сферической каплей. Полная радуга получается вращением угловой функции интенсивности вокруг антисолнечной оси. Для настоящей сцены с дождём нужно суммировать вклад многих капель, но угловая функция одной капли — физическое ядро модели.
