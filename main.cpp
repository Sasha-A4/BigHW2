#include <iostream>
#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <assert.h>
#include <string_view>
#include <optional>
#include <variant>
#include "resurses/fixed.hpp"
#include "resurses/fast_fixed.hpp"
#include "resurses/fixed_oper.hpp"
#include <fstream>
#include <sstream>

using namespace std;

constexpr size_t DEFAULT_ROWS = 36, DEFAULT_COLS = 84;

string readFieldFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cant open file: " + filename);
    }

    string content, line;
    while (getline(file, line)) {
        content += line + '\n';
    }
    return content;
}

#ifndef TYPES
#define TYPES FLOAT,FIXED(31,17),FAST_FIXED(25, 11),FIXED(32, 16),DOUBLE,FAST_FIXED(32, 16)
#endif

#ifndef SIZES
#define SIZES S(36,84), S(100,100)
#endif

struct GridSize {
    size_t rows, cols;
    constexpr GridSize(size_t r, size_t c) : rows(r), cols(c) {}
    constexpr GridSize() : rows(0), cols(0) {}
};

GridSize parseGridSize(const string& sizeStr) {
    size_t start = 2;
    size_t comma = sizeStr.find(',', start);
    size_t end = sizeStr.find(')', comma);
    size_t rows = stoul(sizeStr.substr(start, comma - start));
    size_t cols = stoul(sizeStr.substr(comma + 1, end - comma - 1));

    return GridSize(rows, cols);
}

template<GridSize... Sizes>
struct GridSizesList {
    static constexpr size_t count = sizeof...(Sizes);
    template<size_t I>
    static constexpr GridSize get() {
        constexpr GridSize array[] = {Sizes...};
        return array[I];
    }
};

template<typename List>
bool matchesSize(const GridSize& size) {
    return matchesSizeImpl<List>(size);
}
template<typename List, size_t I = 0>
constexpr bool matchesSizeImpl(const GridSize& size) {
    if constexpr (I >= List::count) {
        return false;
    } else {
        return (List::template get<I>().rows == size.rows && List::template get<I>().cols == size.cols) ||
               matchesSizeImpl<List, I + 1>(size);
    }
}



template<typename NumericType, size_t Rows, size_t Cols>
struct VectorGrid {
    static constexpr std::array<pair<int, int>, 4> directions{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    array<NumericType, 4> vectors[Rows][Cols];

    NumericType &addVector(int x, int y, int dx, int dy, NumericType value) {
        return getVector(x, y, dx, dy) += value;
    }

    NumericType &getVector(int x, int y, int dx, int dy) {
        size_t index = ranges::find(directions, pair(dx, dy)) - directions.begin();
        assert(index < directions.size());
        return vectors[x][y][index];
    }
};

template<
        typename PressureType,
        typename VelocityType,
        typename FlowType,
        size_t Rows = DEFAULT_ROWS,
        size_t Cols = DEFAULT_COLS
>
struct SimulationState {
    PressureType density[256];
    PressureType pressure[Rows][Cols]{};
    PressureType previousPressure[Rows][Cols]{};

    char grid[Rows][Cols + 1];
    int lastUsed[Rows][Cols]{};
    int currentTime{0};
    mt19937 randomGen{1337};
    int directionCount[Rows][Cols]{};

    VectorGrid<VelocityType, Rows, Cols> velocity{};
    VectorGrid<FlowType, Rows, Cols> flowVelocity{};

    SimulationState(const string& gridContent) {
        stringstream ss(gridContent);
        string line;
        size_t row = 0;
        while (getline(ss, line) && row < Rows) {
            size_t col = 0;
            while (col < Cols && col < line.length()) {
                grid[row][col] = line[col];
                col++;
            }
            grid[row][col] = '\0';
            row++;
        }
    }
};

template<typename T>
struct NumericTraits {
    static T fromRaw(int32_t x) { return T(x) / T(1 << 16); }
};
template<size_t N, size_t K>
struct NumericTraits<Fixed<N,K>> {
    static Fixed<N,K> fromRaw(typename Fixed<N,K>::StorageType x) {
        return Fixed<N,K>::fromRaw(x);
    }
};

template<typename T>
struct is_fixed : std::false_type {};
template<size_t N, size_t K>
struct is_fixed<Fixed<N,K>> : std::true_type {};
template<typename T>
inline constexpr bool is_fixed_v = is_fixed<T>::value;
template<typename T>
struct is_fast_fixed : std::false_type {};
template<size_t N, size_t K>
struct is_fast_fixed<FastFixed<N,K>> : std::true_type {};

template<typename T>
inline constexpr bool is_fast_fixed_v = is_fast_fixed<T>::value;

template<typename PressureType, typename VelocityType, typename FlowType, size_t Rows = DEFAULT_ROWS, size_t Cols = DEFAULT_COLS>
class FluidSimulator {
private:
    static constexpr size_t TOTAL_TICKS = 1'000'000;
    static constexpr std::array<pair<int, int>, 4> directions{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

    char grid[Rows][Cols + 1];
    PressureType pressure[Rows][Cols]{}, previousPressure[Rows][Cols];
    int lastUsed[Rows][Cols]{};
    int currentTime = 0;
    mt19937 randomGen{1337};
    PressureType density[256];

    struct VectorGridLocal {
        array<VelocityType, 4> vectors[Rows][Cols]{};

        VelocityType& getVector(int x, int y, int dx, int dy) {
            size_t index = ranges::find(directions, pair(dx, dy)) - directions.begin();
            assert(index < directions.size());
            return vectors[x][y][index];
        }

        VelocityType& addVector(int x, int y, int dx, int dy, VelocityType value) {
            return getVector(x, y, dx, dy) += value;
        }
    };

    VectorGridLocal velocityField{}, flowVelocityField{};
    int directionCount[Rows][Cols]{};

    struct ParticleParameters {
        char type;
        PressureType currentPressure;
        array<VelocityType, 4> vectors;

        void swapWith(FluidSimulator* sim, int x, int y) {
            swap(sim->grid[x][y], type);
            swap(sim->pressure[x][y], currentPressure);
            swap(sim->velocityField.vectors[x][y], vectors);
        }
    };
    template<typename T>
    PressureType convertToPressure(T value) {
        if constexpr (std::is_same_v<T, PressureType>) {
            return value;
        } else {
            return PressureType(value);
        }
    }
    template<typename T>
    VelocityType convertToVelocity(T value) {
        if constexpr (std::is_same_v<T, VelocityType>) {
            return value;
        } else {
            return VelocityType(value);
        }
    }
    template<typename T>
    FlowType convertToFlow(T value) {
        if constexpr (std::is_same_v<T, FlowType>) {
            return value;
        } else {
            return FlowType(value);
        }
    }

    tuple<VelocityType, bool, pair<int, int>> propagateFlow(int x, int y, VelocityType limit) {
        lastUsed[x][y] = currentTime - 1;
        VelocityType totalFlow = 0;
        for (auto [dx, dy] : directions) {
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] != '#' && lastUsed[nx][ny] < currentTime) {
                auto capacity = velocityField.getVector(x, y, dx, dy);
                auto flow = flowVelocityField.getVector(x, y, dx, dy);
                if (flow == capacity) continue;

                auto possibleFlow = min(limit, capacity - flow);
                if (lastUsed[nx][ny] == currentTime - 1) {
                    flowVelocityField.addVector(x, y, dx, dy, possibleFlow);
                    lastUsed[x][y] = currentTime;
                    return {possibleFlow, true, {nx, ny}};
                }
                auto [flowAmount, canPropagate, endPoint] = propagateFlow(nx, ny, possibleFlow);
                totalFlow += flowAmount;
                if (canPropagate) {
                    flowVelocityField.addVector(x, y, dx, dy, flowAmount);
                    lastUsed[x][y] = currentTime;
                    return {flowAmount, canPropagate && endPoint != pair(x, y), endPoint};
                }
            }
        }
        lastUsed[x][y] = currentTime;
        return {totalFlow, false, {0, 0}};
    }

    VelocityType getRandomVelocity() {
        return VelocityType(static_cast<double>(randomGen() & ((1 << 16) - 1)) / (1 << 16));
    }

    void stopPropagation(int x, int y, bool force = false) {
        if (!force) {
            bool shouldStop = true;
            for (auto [dx, dy] : directions) {
                int nx = x + dx, ny = y + dy;
                if (grid[nx][ny] != '#' && lastUsed[nx][ny] < currentTime - 1 && velocityField.getVector(x, y, dx, dy) > VelocityType(0)) {
                    shouldStop = false;
                    break;
                }
            }
            if (!shouldStop) return;
        }
        lastUsed[x][y] = currentTime;
        for (auto [dx, dy] : directions) {
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] == '#' || lastUsed[nx][ny] == currentTime || velocityField.getVector(x, y, dx, dy) > VelocityType(0)) continue;
            stopPropagation(nx, ny);
        }
    }

    VelocityType calculateMoveProbability(int x, int y) {
        VelocityType total = 0;
        for (auto [dx, dy] : directions) {
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] == '#' || lastUsed[nx][ny] == currentTime) continue;
            auto v = velocityField.getVector(x, y, dx, dy);
            if (v < VelocityType(0)) continue;
            total += v;
        }
        return total;
    }

    bool movePropagation(int x, int y, bool isFirst) {
        lastUsed[x][y] = currentTime - isFirst;
        bool moved = false;
        int newX = -1, newY = -1;
        do {
            array<VelocityType, 4> thresholds;
            VelocityType sum = 0;
            for (size_t i = 0; i < directions.size(); ++i) {
                auto [dx, dy] = directions[i];
                int nx = x + dx, ny = y + dy;
                if (grid[nx][ny] == '#' || lastUsed[nx][ny] == currentTime) {
                    thresholds[i] = sum;
                    continue;
                }
                auto v = velocityField.getVector(x, y, dx, dy);
                if (v < VelocityType(0)) {
                    thresholds[i] = sum;
                    continue;
                }
                sum += v;
                thresholds[i] = sum;
            }
            if (sum == VelocityType(0)) break;

            VelocityType randomValue = getRandomVelocity() * sum;
            size_t directionIndex = ranges::upper_bound(thresholds, randomValue) - thresholds.begin();

            auto [dx, dy] = directions[directionIndex];
            newX = x + dx;
            newY = y + dy;
            assert(velocityField.getVector(x, y, dx, dy) > VelocityType(0) && grid[newX][newY] != '#' && lastUsed[newX][newY] < currentTime);

            moved = (lastUsed[newX][newY] == currentTime - 1 || movePropagation(newX, newY, false));
        } while (!moved);

        lastUsed[x][y] = currentTime;
        for (auto [dx, dy] : directions) {
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] != '#' && lastUsed[nx][ny] < currentTime - 1 && velocityField.getVector(x, y, dx, dy) < VelocityType(0)) {
                stopPropagation(nx, ny);
            }
        }
        if (moved && !isFirst) {
            ParticleParameters params{};
            params.swapWith(this, x, y);
            params.swapWith(this, newX, newY);
            params.swapWith(this, x, y);
        }
        return moved;
    }

public:
    FluidSimulator(const SimulationState<PressureType, VelocityType, FlowType, Rows, Cols>& state) {
        memcpy(grid, state.grid, sizeof(grid));
        density[' '] = PressureType(0.01);
        density['.'] = PressureType(1000);

        for (size_t x = 0; x < Rows; ++x) {
            for (size_t y = 0; y < Cols; ++y) {
                if (grid[x][y] == '#') continue;
                for (auto [dx, dy] : directions) {
                    directionCount[x][y] += (grid[x + dx][y + dy] != '#');
                }
            }
        }
    }

    void runSimulation() {
        PressureType gravity = PressureType(0.1);

        for (size_t tick = 0; tick < TOTAL_TICKS; ++tick) {
            PressureType totalPressureChange = 0;

            for (size_t x = 0; x < Rows; ++x) {
                for (size_t y = 0; y < Cols; ++y) {
                    if (grid[x][y] == '#') continue;
                    if (grid[x + 1][y] != '#')
                        velocityField.addVector(x, y, 1, 0, VelocityType(gravity));
                }
            }

            memcpy(previousPressure, pressure, sizeof(pressure));
            for (size_t x = 0; x < Rows; ++x) {
                for (size_t y = 0; y < Cols; ++y) {
                    if (grid[x][y] == '#') continue;
                    for (auto [dx, dy] : directions) {
                        int nx = x + dx, ny = y + dy;
                        if (grid[nx][ny] != '#' && previousPressure[nx][ny] < previousPressure[x][y]) {
                            auto deltaPressure = previousPressure[x][y] - previousPressure[nx][ny];
                            auto force = convertToPressure(deltaPressure);
                            auto& counterForce = velocityField.getVector(nx, ny, -dx, -dy);
                            if (convertToPressure(counterForce) * density[(int)grid[nx][ny]] >= force) {
                                counterForce -= convertToVelocity(force / density[(int)grid[nx][ny]]);
                                continue;
                            }
                            force -= convertToPressure(counterForce) * density[(int)grid[nx][ny]];
                            counterForce = 0;
                            velocityField.addVector(x, y, dx, dy, convertToVelocity(force / density[(int)grid[x][y]]));
                            pressure[x][y] -= force / directionCount[x][y];
                            totalPressureChange -= force / directionCount[x][y];
                        }
                    }
                }
            }

            flowVelocityField = {};
            bool propagationOccurred = false;
            do {
                currentTime += 2;
                propagationOccurred = false;
                for (size_t x = 0; x < Rows; ++x) {
                    for (size_t y = 0; y < Cols; ++y) {
                        if (grid[x][y] != '#' && lastUsed[x][y] != currentTime) {
                            auto [flowAmount, canPropagate, _] = propagateFlow(x, y, VelocityType(1));
                            if (flowAmount > VelocityType(0)) propagationOccurred = true;
                        }
                    }
                }
            } while (propagationOccurred);

            for (size_t x = 0; x < Rows; ++x) {
                for (size_t y = 0; y < Cols; ++y) {
                    if (grid[x][y] == '#') continue;
                    for (auto [dx, dy] : directions) {
                        auto oldVelocity = velocityField.getVector(x, y, dx, dy);
                        auto newVelocity = flowVelocityField.getVector(x, y, dx, dy);
                        if (oldVelocity > VelocityType(0)) {
                            assert(newVelocity <= oldVelocity);
                            velocityField.getVector(x, y, dx, dy) = newVelocity;
                            auto force = convertToPressure(oldVelocity - newVelocity) * density[(int)grid[x][y]];
                            if (grid[x][y] == '.') force *= PressureType(0.8);
                            if (grid[x + dx][y + dy] == '#') {
                                pressure[x][y] += force / directionCount[x][y];
                                totalPressureChange += force / directionCount[x][y];
                            } else {
                                pressure[x + dx][y + dy] += force / directionCount[x + dx][y + dy];
                                totalPressureChange += force / directionCount[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            currentTime += 2;
            propagationOccurred = false;
            for (size_t x = 0; x < Rows; ++x) {
                for (size_t y = 0; y < Cols; ++y) {
                    if (grid[x][y] != '#' && lastUsed[x][y] != currentTime) {
                        if (getRandomVelocity() < calculateMoveProbability(x, y)) {
                            propagationOccurred = true;
                            movePropagation(x, y, true);
                        } else {
                            stopPropagation(x, y, true);
                        }
                    }
                }
            }

            if (propagationOccurred) {
                cout << "Tick " << tick << ":\n";
                for (size_t x = 0; x < Rows; ++x) {
                    cout << grid[x] << "\n";
                }
            }
        }
    }
};

template<typename P, typename V, typename VF>
void executeSimulation(size_t rows, size_t cols, const string& gridData) {
    SimulationState<P, V, VF, DEFAULT_ROWS, DEFAULT_COLS> state(gridData);
    FluidSimulator<P, V, VF, DEFAULT_ROWS, DEFAULT_COLS> simulator(state);
    simulator.runSimulation();
}

template<typename T>
std::string getTypeName() {
    if constexpr (std::is_same_v<T, float>) {
        return "float";
    } else if constexpr (std::is_same_v<T, double>) {
        return "double";
    } else if constexpr (is_fixed_v<T>) {
        return "Fixed<" + std::to_string(T::bits) + "," + std::to_string(T::frac) + ">";
    } else if constexpr (is_fast_fixed_v<T>) {
        return "FastFixed<" + std::to_string(T::bits) + "," + std::to_string(T::frac) + ">";
    } else {
        return "unknown";
    }
}

pair<size_t, size_t> parseFixedParameters(const string& typeStr) {
    size_t start = typeStr.find('(') + 1;
    size_t comma = typeStr.find(',', start);
    size_t end = typeStr.find(')', comma);

    size_t bits = stoul(typeStr.substr(start, comma - start));
    size_t frac = stoul(typeStr.substr(comma + 1, end - comma - 1));
    return {bits, frac};
}

template<typename T>
static bool isMatchingType(const string& typeStr) {
    if constexpr (std::is_same_v<T, float>) {
        return typeStr == "FLOAT";
    } else if constexpr (std::is_same_v<T, double>) {
        return typeStr == "DOUBLE";
    } else if constexpr (is_fixed_v<T>) {
        if (!typeStr.starts_with("FIXED(")) return false;
        auto [bits, frac] = parseFixedParameters(typeStr);
        return bits == T::bits && frac == T::frac;
    } else if constexpr (is_fast_fixed_v<T>) {
        if (!typeStr.starts_with("FAST_FIXED(")) return false;
        auto [bits, frac] = parseFixedParameters(typeStr);
        return bits == T::bits && frac == T::frac;
    }
    return false;
}
template<typename... Types>
struct TypeList {
    static constexpr size_t count = sizeof...(Types);
    template<size_t I>
    using type_at = typename std::tuple_element<I, std::tuple<Types...>>::type;
};

template<typename AllowedTypes, typename SelectedTypes>
struct TypeSelector {
    template<typename... Selected>
    static bool tryCombinations(const string& pType, const string& vType, const string& flowType,
                                const GridSize& size, const string& gridData) {
        return tryPressureTypes<0>(pType, vType, flowType, size, gridData);
    }

private:
    template<size_t I>
    static bool tryPressureTypes(const string& pType, const string& vType, const string& flowType,
                                 const GridSize& size, const string& gridData) {
        if constexpr (I >= AllowedTypes::count) {
            return false;
        } else {
            using Pressure = typename AllowedTypes::template type_at<I>;
            return tryWithPressureType<Pressure>(pType, vType, flowType, size, gridData) ||
                   tryPressureTypes<I + 1>(pType, vType, flowType, size, gridData);
        }
    }
    template<typename Pressure>
    static bool tryWithPressureType(const string& pType, const string& vType, const string& flowType,
                                    const GridSize& size, const string& gridData) {
        if (!isMatchingType<Pressure>(pType)) return false;
        return tryVelocityTypes<Pressure, 0>(pType, vType, flowType, size, gridData);
    }
    template<typename Pressure, size_t I>
    static bool tryVelocityTypes(const string& pType, const string& vType, const string& flowType,
                                 const GridSize& size, const string& gridData) {
        if constexpr (I >= AllowedTypes::count) {
            return false;
        } else {
            using Velocity = typename AllowedTypes::template type_at<I>;
            return tryWithVelocityType<Pressure, Velocity>(pType, vType, flowType, size, gridData) ||
                   tryVelocityTypes<Pressure, I + 1>(pType, vType, flowType, size, gridData);
        }
    }
    template<typename Pressure, typename Velocity>
    static bool tryWithVelocityType(const string& pType, const string& vType, const string& flowType,
                                    const GridSize& size, const string& gridData) {
        if (!isMatchingType<Velocity>(vType)) return false;
        return tryFlowTypes<Pressure, Velocity, 0>(pType, vType, flowType, size, gridData);
    }
    template<typename Pressure, typename Velocity, size_t I>
    static bool tryFlowTypes(const string& pType, const string& vType, const string& flowType,
                             const GridSize& size, const string& gridData) {
        if constexpr (I >= AllowedTypes::count) {
            return false;
        } else {
            using Flow = typename AllowedTypes::template type_at<I>;
            return tryWithFlowType<Pressure, Velocity, Flow>(pType, vType, flowType, size, gridData) ||
                   tryFlowTypes<Pressure, Velocity, I + 1>(pType, vType, flowType, size, gridData);
        }
    }
    template<typename Pressure, typename Velocity, typename Flow>
    static bool tryWithFlowType(const string& pType, const string& vType, const string& flowType,
                                const GridSize& size, const string& gridData) {
        if (!isMatchingType<Flow>(flowType)) return false;



        executeSimulation<Pressure, Velocity, Flow>(size.rows, size.cols, gridData);
        return true;
    }
};

template<typename... Types>
bool tryAllTypeCombinations(const string& pType, const string& vType, const string& flowType,
                            const GridSize& size, const string& gridData) {
    return TypeSelector<TypeList<Types...>, TypeList<>>::tryCombinations(pType, vType, flowType, size, gridData);
}

bool createAndRunSimulation(const string& pType, const string& vType, const string& flowType,
                            const GridSize& size, const string& gridData) {
    try {


#define S(r, c) GridSize(r, c)
        using GridSizesListType = GridSizesList<SIZES>;
#undef S

        if (!matchesSize<GridSizesListType>(size)) {
            cout << "size error" << endl;
            return false;
        }

#define FLOAT float
#define DOUBLE double
#define FIXED(r, c) Fixed<r, c>
#define FAST_FIXED(r, c) FastFixed<r, c>
        if (!tryAllTypeCombinations<TYPES>(pType, vType, flowType, size, gridData)) {
            cout << "type combination error" << endl;
            return false;
        }
        return true;
    }
    catch (const exception& e) {
        cout << "error: " << e.what() << endl;
        return false;
    }
}

string getArgument(string_view argName, int argc, char** argv, string_view defaultValue) {
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == argName) {
            return argv[i + 1];
        }
    }
    return string(defaultValue);
}

bool isValidType(const string& typeStr) {
    if (typeStr == "FLOAT" || typeStr == "DOUBLE") return true;

    if (typeStr.starts_with("FIXED(") || typeStr.starts_with("FAST_FIXED(")) {
        size_t commaPos = typeStr.find(',');
        if (commaPos == string::npos) return false;

        size_t endPos = typeStr.find(')', commaPos);
        if (endPos == string::npos) return false;

        return true;
    }

    return false;
}

GridSize determineFieldSize(const string& gridData) {
    size_t rows = 0, cols = 0;
    stringstream ss(gridData);
    string line;

    while (getline(ss, line)) {
        if (cols == 0) cols = line.length();
        else if (line.length() != cols) {
            throw runtime_error("field format is wrong");
        }
        rows++;
    }


    return GridSize(rows, cols);
}

int main(int argc, char** argv) {
    string pressureType = getArgument("--p-type", argc, argv, "FAST_FIXED(32, 16)");
    string velocityType = getArgument("--v-type", argc, argv, "FIXED(32, 16)");
    string flowType = getArgument("--v-flow-type", argc, argv, "FAST_FIXED(32, 16)");
    string fieldFile = getArgument("--field", argc, argv, "../main.txt");
    string gridData;
    gridData = readFieldFromFile(fieldFile);
    GridSize size;
    size = determineFieldSize(gridData);


    if (!isValidType(pressureType)) {
        cout << "Invalid p_type: " << pressureType << endl;
        return 1;
    }
    if (!isValidType(velocityType)) {
        cout << "Invalid v_type: " << velocityType << endl;
        return 1;
    }
    if (!isValidType(flowType)) {
        cout << "Invalid v_flow_type: " << flowType << endl;
        return 1;
    }
    if (!createAndRunSimulation(pressureType, velocityType, flowType, size, gridData)) {
        cout << "simulation error" << endl;
        return 1;
    }

    return 0;
}

