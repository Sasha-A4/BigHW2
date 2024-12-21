#include <iostream>
#include <array>
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <fstream>
#include <filesystem>

using namespace std;

constexpr size_t ROWS = 36, COLS = 84;
constexpr size_t STEPS = 1'000'000;
constexpr std::array<pair<int, int>, 4> DIRECTIONS{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};


char grid[ROWS][COLS + 1];

void loadGridFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "can't open " << filename << endl;
        exit(1);
    }

    string line;
    for (size_t i = 0; i < ROWS; ++i) {
        getline(file, line);
        strcpy(grid[i], line.c_str());
    }

    file.close();
}

struct FixedPoint {
    constexpr FixedPoint(int val): value(val << 16) {}
    constexpr FixedPoint(float f): value(f * (1 << 16)) {}
    constexpr FixedPoint(double f): value(f * (1 << 16)) {}
    constexpr FixedPoint(): value(0) {}

    static constexpr FixedPoint from_raw(int32_t x) {
        FixedPoint result;
        result.value = x;
        return result;
    }

    int32_t value;

    auto operator<=>(const FixedPoint&) const = default;
    bool operator==(const FixedPoint&) const = default;
};

//static constexpr FixedPoint INF = FixedPoint::from_raw(std::numeric_limits<int32_t>::max());
static constexpr FixedPoint EPSILON = FixedPoint::from_raw(DIRECTIONS.size());

FixedPoint operator+(FixedPoint a, FixedPoint b) {
    return FixedPoint::from_raw(a.value + b.value);
}

FixedPoint operator-(FixedPoint a, FixedPoint b) {
    return FixedPoint::from_raw(a.value - b.value);
}

FixedPoint operator*(FixedPoint a, FixedPoint b) {
    return FixedPoint::from_raw(((int64_t) a.value * b.value) >> 16);
}

FixedPoint operator/(FixedPoint a, FixedPoint b) {
    return FixedPoint::from_raw(((int64_t) a.value << 16) / b.value);
}

FixedPoint &operator+=(FixedPoint &a, FixedPoint b) {
    return a = a + b;
}

FixedPoint &operator-=(FixedPoint &a, FixedPoint b) {
    return a = a - b;
}

FixedPoint &operator*=(FixedPoint &a, FixedPoint b) {
    return a = a * b;
}

FixedPoint &operator/=(FixedPoint &a, FixedPoint b) {
    return a = a / b;
}

FixedPoint operator-(FixedPoint x) {
    return FixedPoint::from_raw(-x.value);
}

FixedPoint abs(FixedPoint x) {
    if (x.value < 0) {
        x.value = -x.value;
    }
    return x;
}

ostream &operator<<(ostream &out, FixedPoint x) {
    return out << x.value / (double) (1 << 16);
}

FixedPoint rhoTable[256];

FixedPoint pressure[ROWS][COLS]{}, previousPressure[ROWS][COLS];

struct VectorField {
    array<FixedPoint, DIRECTIONS.size()> velocity[ROWS][COLS];
    FixedPoint &add(int x, int y, int dx, int dy, FixedPoint delta) {
        return get(x, y, dx, dy) += delta;
    }

    FixedPoint &get(int x, int y, int dx, int dy) {
        size_t idx = ranges::find(DIRECTIONS, pair(dx, dy)) - DIRECTIONS.begin();
        assert(idx < DIRECTIONS.size());
        return velocity[x][y][idx];
    }
};

VectorField velocityField{}, flowField{};
int lastUsage[ROWS][COLS]{};
int updateTick = 0;

mt19937 randomGen(1337);

tuple<FixedPoint, bool, pair<int, int>> propagateFlow(int x, int y, FixedPoint limit) {
    lastUsage[x][y] = updateTick - 1;
    FixedPoint result = 0;
    for (auto [dx, dy] : DIRECTIONS) {
        int nx = x + dx, ny = y + dy;
        if (grid[nx][ny] != '#' && lastUsage[nx][ny] < updateTick) {
            auto capacity = velocityField.get(x, y, dx, dy);
            auto flow = flowField.get(x, y, dx, dy);
            if (flow == capacity) {
                continue;
            }
            auto delta = min(limit, capacity - flow);
            if (lastUsage[nx][ny] == updateTick - 1) {
                flowField.add(x, y, dx, dy, delta);
                lastUsage[x][y] = updateTick;
                return {delta, 1, {nx, ny}};
            }
            auto [t, propagate, end] = propagateFlow(nx, ny, delta);
            result += t;
            if (propagate) {
                flowField.add(x, y, dx, dy, t);
                lastUsage[x][y] = updateTick;
                return {t, propagate && end != pair(x, y), end};
            }
        }
    }
    lastUsage[x][y] = updateTick;
    return {result, 0, {0, 0}};
}

FixedPoint generateRandom() {
    return FixedPoint::from_raw((randomGen() & ((1 << 16) - 1)));
}

void propagateStop(int x, int y, bool force = false) {
    if (!force) {
        bool stop = true;
        for (auto [dx, dy] : DIRECTIONS) {
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] != '#' && lastUsage[nx][ny] < updateTick - 1 && velocityField.get(x, y, dx, dy) > 0) {
                stop = false;
                break;
            }
        }
        if (!stop) {
            return;
        }
    }
    lastUsage[x][y] = updateTick;
    for (auto [dx, dy] : DIRECTIONS) {
        int nx = x + dx, ny = y + dy;
        if (grid[nx][ny] == '#' || lastUsage[nx][ny] == updateTick || velocityField.get(x, y, dx, dy) > 0) {
            continue;
        }
        propagateStop(nx, ny);
    }
}

FixedPoint calculateMoveProbability(int x, int y) {
    FixedPoint sum = 0;
    for (size_t i = 0; i < DIRECTIONS.size(); ++i) {
        auto [dx, dy] = DIRECTIONS[i];
        int nx = x + dx, ny = y + dy;
        if (grid[nx][ny] == '#' || lastUsage[nx][ny] == updateTick) {
            continue;
        }
        auto v = velocityField.get(x, y, dx, dy);
        if (v < 0) {
            continue;
        }
        sum += v;
    }
    return sum;
}

struct Particle {
    char type;
    FixedPoint currentPressure;
    array<FixedPoint, DIRECTIONS.size()> velocity;

    void swapWith(int x, int y) {
        swap(grid[x][y], type);
        swap(pressure[x][y], currentPressure);
        swap(velocityField.velocity[x][y], velocity);
    }
};

bool propagateMove(int x, int y, bool isFirstMove) {
    lastUsage[x][y] = updateTick - isFirstMove;
    bool result = false;
    int nx = -1, ny = -1;
    do {
        std::array<FixedPoint, DIRECTIONS.size()> temp;
        FixedPoint sum = 0;
        for (size_t i = 0; i < DIRECTIONS.size(); ++i) {
            auto [dx, dy] = DIRECTIONS[i];
            int nx = x + dx, ny = y + dy;
            if (grid[nx][ny] == '#' || lastUsage[nx][ny] == updateTick) {
                temp[i] = sum;
                continue;
            }
            auto v = velocityField.get(x, y, dx, dy);
            if (v < 0) {
                temp[i] = sum;
                continue;
            }
            sum += v;
            temp[i] = sum;
        }

        if (sum == 0) {
            break;
        }

        FixedPoint randomPressure = generateRandom() * sum;
        size_t dirIdx = std::ranges::upper_bound(temp, randomPressure) - temp.begin();

        auto [dx, dy] = DIRECTIONS[dirIdx];
        nx = x + dx;
        ny = y + dy;
        assert(velocityField.get(x, y, dx, dy) > 0 && grid[nx][ny] != '#' && lastUsage[nx][ny] < updateTick);

        result = (lastUsage[nx][ny] == updateTick - 1 || propagateMove(nx, ny, false));
    } while (!result);
    lastUsage[x][y] = updateTick;
    for (size_t i = 0; i < DIRECTIONS.size(); ++i) {
        auto [dx, dy] = DIRECTIONS[i];
        int nx = x + dx, ny = y + dy;
        if (grid[nx][ny] != '#' && lastUsage[nx][ny] < updateTick - 1 && velocityField.get(x, y, dx, dy) < 0) {
            propagateStop(nx, ny);
        }
    }
    if (result) {
        if (!isFirstMove) {
            Particle particle{};
            particle.swapWith(x, y);
            particle.swapWith(nx, ny);
            particle.swapWith(x, y);
        }
    }
    return result;
}

int directions[ROWS][COLS]{};

int main() {
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
    freopen("main.txt", "r", stdin);
    loadGridFromFile("main.txt");

    rhoTable[' '] = 0.01;
    rhoTable['.'] = 1000;
    FixedPoint gravity = 0.1;

    for (size_t x = 0; x < ROWS; ++x) {
        for (size_t y = 0; y < COLS; ++y) {
            if (grid[x][y] == '#')
                continue;
            for (auto [dx, dy] : DIRECTIONS) {
                directions[x][y] += (grid[x + dx][y + dy] != '#');
            }
        }
    }

    for (size_t i = 0; i < STEPS; ++i) {

        FixedPoint totalPressureChange = 0;
        for (size_t x = 0; x < ROWS; ++x) {
            for (size_t y = 0; y < COLS; ++y) {
                if (grid[x][y] == '#')
                    continue;
                if (grid[x + 1][y] != '#')
                    velocityField.add(x, y, 1, 0, gravity);
            }
        }

        memcpy(previousPressure, pressure, sizeof(pressure));
        for (size_t x = 0; x < ROWS; ++x) {
            for (size_t y = 0; y < COLS; ++y) {
                if (grid[x][y] == '#')
                    continue;
                for (auto [dx, dy] : DIRECTIONS) {
                    int nx = x + dx, ny = y + dy;
                    if (grid[nx][ny] != '#' && previousPressure[nx][ny] < previousPressure[x][y]) {
                        auto deltaPressure = previousPressure[x][y] - previousPressure[nx][ny];
                        auto force = deltaPressure;
                        auto &contr = velocityField.get(nx, ny, -dx, -dy);
                        if (contr * rhoTable[(int) grid[nx][ny]] >= force) {
                            contr -= force / rhoTable[(int) grid[nx][ny]];
                            continue;
                        }
                        force -= contr * rhoTable[(int) grid[nx][ny]];
                        contr = 0;
                        velocityField.add(x, y, dx, dy, force / rhoTable[(int) grid[x][y]]);
                        pressure[x][y] -= force / directions[x][y];
                        totalPressureChange -= force / directions[x][y];
                    }
                }
            }
        }

        flowField = {};
        bool propagated = false;
        do {
            updateTick += 2;
            propagated = 0;
            for (size_t x = 0; x < ROWS; ++x) {
                for (size_t y = 0; y < COLS; ++y) {
                    if (grid[x][y] != '#' && lastUsage[x][y] != updateTick) {
                        auto [t, localProp, _] = propagateFlow(x, y, 1);
                        if (t > 0) {
                            propagated = 1;
                        }
                    }
                }
            }
        } while (propagated);

        for (size_t x = 0; x < ROWS; ++x) {
            for (size_t y = 0; y < COLS; ++y) {
                if (grid[x][y] == '#')
                    continue;
                for (auto [dx, dy] : DIRECTIONS) {
                    auto oldVelocity = velocityField.get(x, y, dx, dy);
                    auto newVelocity = flowField.get(x, y, dx, dy);
                    if (oldVelocity > 0) {
                        assert(newVelocity <= oldVelocity);
                        velocityField.get(x, y, dx, dy) = newVelocity;
                        auto force = (oldVelocity - newVelocity) * rhoTable[(int) grid[x][y]];
                        if (grid[x][y] == '.')
                            force *= 0.8;
                        if (grid[x + dx][y + dy] == '#') {
                            pressure[x][y] += force / directions[x][y];
                            totalPressureChange += force / directions[x][y];
                        } else {
                            pressure[x + dx][y + dy] += force / directions[x + dx][y + dy];
                            totalPressureChange += force / directions[x + dx][y + dy];
                        }
                    }
                }
            }
        }

        updateTick += 2;
        propagated = false;
        for (size_t x = 0; x < ROWS; ++x) {
            for (size_t y = 0; y < COLS; ++y) {
                if (grid[x][y] != '#' && lastUsage[x][y] != updateTick) {
                    if (generateRandom() < calculateMoveProbability(x, y)) {
                        propagated = true;
                        propagateMove(x, y, true);
                    } else {
                        propagateStop(x, y, true);
                    }
                }
            }
        }

        if (propagated) {
            cout << "Tick number " << i << ":\n";
            for (size_t x = 0; x < ROWS; ++x) {
                cout << grid[x] << "\n";
            }
        }
    }
}
