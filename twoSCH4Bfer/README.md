# twoSCH4BferLib — OpenFOAM 14 (openfoam.org) библиотека горения CH₄/O₂

Библиотека добавляет в OpenFOAM 14 два компонента:

1. **Транспортный пакет `twoSCH4Bfer`** — степенная вязкость с температурной
   зависимостью и постоянным числом Прандтля (или постоянной теплопроводностью).
2. **Скорость реакции `ArrheniusReactionRateExtend`** — закон Аррениуса с
   дополнительным множителем $f(\varphi)$, зависящим от коэффициента
   эквивалентности $\varphi$ (двухстадийный механизм горения метана).

> Это порт самописной библиотеки из ESI/OpenCFD (`com`) 2506 в
> openfoam.org 14. Полный набор уравнений состояния (19 комбинаций) и полная
> регистрация chemistry/reduction-методов уже включены.

---

## 1. Предварительные требования

- **OpenFOAM 14 (openfoam.org)**. Проверка:
  ```bash
  source /opt/openfoam14/etc/bashrc
  echo $WM_PROJECT_VERSION   # → 14
  ```
- Компилятор `g++` (стандарт C++14, задаётся автоматически `wmake`).
- Права на запись в `$FOAM_USER_LIBBIN`
  (по умолчанию `~/OpenFOAM/user-14/platforms/linux64GccDPInt32Opt/lib`).

> ⚠️ **Важно про окружение.** В `~/.bashrc` этого компьютера загружен
> **OpenFOAM v2506** (`source ~/OpenFOAM/OpenFOAM-v2506/etc/bashrc`). Поэтому
> обычный `wmake` будет собирать **не ту** версию. Всегда собирайте в чистом
> окружении (см. раздел 3).

---

## 2. Состав библиотеки

```
src/twoSCH4Bfer/
├── twoSCH4Bfer.H / .C / .I.H          # транспортный пакет (mu(T)=mu0*(T/T0)^alpha)
├── ArrheniusReactionRateExtend.H/.I.H # скорость реакции + множитель f(phi)
├── forTwoSCH4Bfer.H                   # макрокаскады (все EoS/thermo комбинации)
├── twoSCH4BferThermos.C               # регистрация thermo (psi + rhoFluid)
├── makeReactionsTwoS.C                # регистрация реакций
├── makeChemistryModelsTwoS.C          # регистрация Standard chemistry model
├── makeChemistryReductionMethodsTwoS.C# регистрация reduction-методов
└── Make/
    ├── files
    └── options
```

Собранный результат — одна динамическая библиотека:

```
$FOAM_USER_LIBBIN/twoSCH4BferLib.so
```

---

## 3. Компиляция (установка)

```bash
# 1) Чистое окружение + org 14 (важно из-за v2506 в ~/.bashrc!)
env -i HOME="$HOME" USER="$USER" TERM="$TERM" \
  PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
  bash --noprofile --norc -c \
  'source /opt/openfoam14/etc/bashrc && \
   cd "$HOME/OpenFOAM/OpenFOAM Driver/OpenFOAM_org/src/twoSCH4Bfer" && \
   wmake libso'
```

Что делает `wmake libso`:
- `wmakeLnInclude` — генерирует `lnInclude/` (символические ссылки на `.H`);
- компилирует 4 `.C` файла в `Make/linux64GccDPInt32Opt/*.o`;
- линкует `twoSCH4BferLib.so` в `$FOAM_USER_LIBBIN`.

Проверка, что библиотека собралась:

```bash
ls -la "$HOME/OpenFOAM/user-t6-14/platforms/linux64GccDPInt32Opt/lib/twoSCH4BferLib.so"
```

### Пересборка после правок

```bash
# то же самое окружение, затем:
wmake libso
```

`wmake` сам перекомпилирует только изменённые исходники. Для **полной чистой**
пересборки:

```bash
wclean && wmake libso
```

---

## 4. Подключение библиотеки к кейсу

В `system/controlDict` кейса:

```cpp
libs
(
    "$FOAM_USER_LIBBIN/twoSCH4BferLib.so"
);
```

> Если используете кейс в другом месте — `$FOAM_USER_LIBBIN` должно
> указывать на ту версию OpenFOAM, которой вы запускаете solver. При запуске
> с org 14 это `~/OpenFOAM/user-t6-14/...`.

---

## 5. Настройка thermo (constant/physicalProperties)

```cpp
thermoType
{
    type            hePsiThermo;      // psi-based (p, T)
    mixture         multicomponentMixture;
    transport       twoSCH4Bfer;
    thermo          janaf;            // или hConst / eConst
    energy          sensibleEnthalpy; // или sensibleInternalEnergy
    equationOfState perfectGas;       // см. список ниже
    specie          specie;
}

defaultSpecie N2;   // вид, который рассчитывается как 1 - сумма остальных

#include "thermo.compressibleGas"
```

> ⚠️ **Ключевой момент org 14:** в словаре пишется `mixture multicomponentMixture;`
> (а **не** `coefficientMulticomponentMixture`). Причина: класс
> `coefficientMulticomponentMixture` в org 14 имеет `typeName()` =
> `multicomponentMixture<...>`. Это совпадает с тем, что ожидает RTS-таблица
> `fluidMulticomponentThermo`.

### Транспортные коэффициенты (в `thermo.compressibleGas`)

Каждому виду задаётся блок `transport`:

```cpp
O2
{
    specie      { molWeight 31.9988; }
    thermodynamics { ... janaf полиномы ... }
    transport
    {
        mu          1.8405e-5;   // [Pa.s] опорная вязкость mu0
        T0          300;         // [K]    опорная температура
        muExponent  0.6759;      // [-]    показатель alpha: mu = mu0*(T/T0)^alpha
        Pr          0.7;         // [-]    ИЛИ kappa <W/m/K> (одно из двух!)
    }
}
```

Правило выбора `Pr` или `kappa`:
- задан `Pr`  → теплопроводность считается как $\kappa = C_p \mu / Pr$;
- задан `kappa` → теплопроводность постоянна;
- **одновременно** `Pr` и `kappa` задавать нельзя (FatalIOError).

### Доступные уравнения состояния (полный набор)

| equationOfState | mixture (multicomponent) | single-component |
|---|---|---|
| `perfectGas` | ✅ | ✅ |
| `PengRobinsonGas` | ✅ | ✅ |
| `incompressiblePerfectGas` | ✅ | ✅ |
| `rPolynomial` | ✅ | ✅ |
| `rhoConst` | ✅ | ✅ |
| `perfectFluid` | ❌ (noCoefficientMixing) | ✅ |
| `adiabaticPerfectFluid` | ❌ (noCoefficientMixing) | ✅ |
| `icoPolynomial` (+ `hPolynomialThermo`) | ❌ (noCoefficientMixing) | ✅ |

> Причина ограничения: org 14 помечает `perfectFluid`, `adiabaticPerfectFluid`,
> `icoPolynomial` как `noCoefficientMixing` — их нельзя использовать в
> многокомпонентной смеси. Для горения CH₄/O₂/N₂ достаточно `perfectGas`.

---

## 6. Реакции (constant/reactions)

Двухстадийный механизм с множителем $f(\varphi)$:

```cpp
reactions
{
    R1
    {
        type        irreversibleArrheniusTwoSCH4Bfer;
        reaction    "CH4^0.5 + 1.5O2^0.65 = CO + 2H2O";
        A           1.8677e9;
        beta        0;
        Ta          17866;
    }
    R2
    {
        type        reversibleArrheniusTwoSCH4Bfer;
        reaction    "CO^1 + 0.5O2^0.5 = CO2";
        A           1.408e6;
        beta        0.8;
        Ta          6039;
    }
}
```

### Доступные типы реакций

| тип | назначение |
|---|---|
| `irreversibleArrheniusTwoSCH4Bfer` | необратимая |
| `reversibleArrheniusTwoSCH4Bfer` | обратимая |
| `nonEquilibriumReversibleArrheniusTwoSCH4Bfer` | не-равновесная обратимая |

### Откуда берётся $\varphi$

Множитель $f(\varphi)$ читает коэффициент эквивалентности один раз за запуск из
файла **`constant/phiEq`**:

```cpp
FoamFile { ... class dictionary; object phiEq; }
value   0.8;
```

Если файл отсутствует или нет записи `value` — FatalError при старте.

> Внутренняя логика: по имени реакции (`R1` или `R2`) выбирается своя функция
> $f_1(\varphi)$ или $f_2(\varphi)$ (см. `ArrheniusReactionRateExtendI.H`).
> Реакции должны называться именно `R1` / `R2`, иначе множитель $f=1$.

---

## 7. Прочие файлы кейса (org 14)

| файл | назначение |
|---|---|
| `constant/chemistryProperties` | `initialChemicalTimeStep` + блок `ode` (solver, tol) |
| `constant/reactionProperties` | `reactionModel laminar;` |
| `constant/momentumTransport` | `simulationType laminar;` |
| `constant/thermophysicalTransport` | **число Льюиса** (см. ниже) |
| `system/controlDict` | `solver multicomponentFluid;` + `libs` |

### Число Льюиса Le = 1

В org 14 ламинарная модель тепломассопереноса **по умолчанию** —
`unityLewisFourier` (Le = 1), т.е. файл `constant/thermophysicalTransport`
можно не создавать. Чтобы задать явно:

```cpp
laminar
{
    model           unityLewisFourier;   // Le = 1: D = kappa/Cp = alpha
}
```

Альтернативные ламинарные модели (если понадобится Le ≠ 1):
`FickianFourier`, `Fickian`, `MaxwellStefanFourier`, `Fourier`.

---

## 8. Запуск solver

org 14 использует **модульный** solver `foamRun` (модуль `multicomponentFluid`):

```bash
# то же чистое окружение org 14
env -i HOME="$HOME" USER="$USER" TERM="$TERM" \
  PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
  bash --noprofile --norc -c \
  'source /opt/openfoam14/etc/bashrc && \
   cd <путь_к_кейсу> && foamRun'
```

**Важно: путь к кейсу не должен содержать пробелов.** OpenFOAM (org 14) падает
на `fileName::stripInvalid()` при наличии пробела в пути. Если кейс лежит в
`.../OpenFOAM Driver/...`, скопируйте его в путь без пробелов (например
`~/OpenFOAM/run/<case>`).

---

## 9. Различия com 2506 → org 14 (важно при переносе кейсов)

| аспект | com 2506 | org 14 (этот порт) |
|---|---|---|
| solver | `reactingFoam` / `reactingFoamUnityLe` | `foamRun` + модуль `multicomponentFluid` |
| thermo файл | `constant/thermophysicalProperties` | `constant/physicalProperties` |
| mixture | `reactingMixture` | `multicomponentMixture` |
| chemistry reader | `foamChemistryReader` + `foamChemistryFile` | не нужен (всё через `#include "reactions"`) |
| combustion | `constant/combustionProperties` | `constant/reactionProperties` (`reactionModel laminar`) |
| turbulence | `constant/turbulenceProperties` | `constant/momentumTransport` |
| тип реакции | `...ArrheniusExtendReaction` | `...ArrheniusTwoSCH4Bfer` (без суффикса `Reaction`) |
| TDAC | отдельный `TDACChemistryModel` | объединён в `Standard` chemistry model |
| reduction/tabulation | отдельные `makeChemistryReduction...` | включены в эту либу (`makeChemistryReductionMethodsTwoS.C`) |

---

## 10. Частые ошибки и их решение

### 10.1 `Unknown mixture type multicomponentMixture`
Библиотека не загрузилась. Проверьте:
- в `controlDict` подключена `libs (... twoSCH4BferLib.so ...)`;
- библиотека собрана именно для org 14 (не v2506).

### 10.2 `undefined symbol: chemistryReductionMethod<...>::dictionaryConstructorTablePtr_`
Старая сборка без `makeChemistryReductionMethodsTwoS.C`. Пересоберите либу
(файл уже добавлен в `Make/files`).

### 10.3 `Min sum of specie mass fractions ... < 0.999`
В org 14 сумма массовых долей считается **со 2-го вида**, первый вид в списке
`species` должен иметь начальную долю ≈ 0 (обычно это продукт, например `CO`),
а сумма остальных (реагенты + разбавитель) = 1.

### 10.4 `fileName::stripInvalid() called for invalid fileName`
В пути кейса есть пробел (или недопустимый символ). Перенесите кейс в путь без
пробелов.

### 10.5 `Cannot open .../constant/phiEq`
Создайте `constant/phiEq` с записью `value <число>;`.

---

## 11. Быстрый старт (минимальный кейс)

```bash
# 1) собрать либу (раздел 3)
# 2) подготовить кейс:
#    - constant/physicalProperties (thermoType + transport twoSCH4Bfer)
#    - constant/thermo.compressibleGas (виды + transport {mu T0 muExponent Pr})
#    - constant/reactions (R1/R2 с типом ...ArrheniusTwoSCH4Bfer)
#    - constant/phiEq (value 0.8)
#    - constant/chemistryProperties / reactionProperties / momentumTransport
#    - system/controlDict (solver multicomponentFluid + libs)
#    - 0/ поля: p, U, T, CH4, O2, N2, CO, CO2, H2O, Ydefault
# 3) запустить: foamRun (в чистом окружении org 14)
```

---

## 12. Лицензия / происхождение

Портировано из ESI/OpenCFD (com) 2506 в openfoam.org 14.
Механизм реакций основан на двухстадийной модели Bui-Pham (1992), модифицирован
для зависимости скорости от коэффициента эквивалентности $\varphi$.
