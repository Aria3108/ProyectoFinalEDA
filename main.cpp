#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <unordered_map>
using namespace std;

struct Node {
    long long id;
    long double lat;
    long double lon;
};

struct Edge {
    long long id;
    Node* start;
    Node* end;
};

struct Point {
    long double x, y;
    Point() : x(0), y(0) {}
    Point(long double x_, long double y_) : x(x_), y(y_) {}
    bool operator==(Point other) const { return x == other.x && y == other.y; }
};

struct MBR {
    Point min, max;
    MBR() : min(), max() {}
    MBR(Point min_, Point max_) : min(min_), max(max_) {}
};

bool menorx(const Point& a, const Point& b) { return a.x < b.x; }
bool menory(const Point& a, const Point& b) { return a.y < b.y; }

bool puntoDentroDeMBR(Point p, MBR mbr) {
    return (p.x >= mbr.min.x && p.x <= mbr.max.x &&
            p.y >= mbr.min.y && p.y <= mbr.max.y);
}

// cambiar || por && si explota
bool Contenida(Edge ar, MBR mbr) {
    if (puntoDentroDeMBR({ ar.start->lon, ar.start->lat }, mbr) && puntoDentroDeMBR({ ar.end->lon, ar.end->lat }, mbr)) {
        return true;
    }

    Point p1 = mbr.min, p2 = mbr.max;
    Point pa1 = { ar.start->lon, ar.start->lat }, pa2 = { ar.end->lon, ar.end->lat };
    long double A = pa2.y - pa1.y;
    long double B = -(pa2.x - pa1.x);
    long double C = pa2.x * pa1.y - pa1.x * pa2.y;
    vector<Point> intersections;

    if (B != 0) {
        long double y = (-C - A * p1.x) / B;
        if (p1.y <= y && y <= p2.y) {
            intersections.push_back({ p1.x, y });
        }
    }

    if (B != 0) {
        long double y = (-C - A * p2.x) / B;
        if (p1.y <= y && y <= p2.y) {
            intersections.push_back({ p2.x, y });
        }
    }

    if (A != 0) {
        long double x = (-C - B * p1.y) / A;
        if (p1.x <= x && x <= p2.x) {
            intersections.push_back({ x, p1.y });
        }
    }

    if (A != 0) {
        long double x = (-C - B * p2.y) / A;
        if (p1.x <= x && x <= p2.x) {
            intersections.push_back({ x, p2.y });
        }
    }

    if (intersections.size() == 2) {
        Point i1 = intersections[0];
        Point i2 = intersections[1];
        if (i1 == i2) {
            return false;
        } else if (puntoDentroDeMBR(i1, mbr) && puntoDentroDeMBR(i2, mbr)) {
            if ((i1.x >= min(pa1.x, pa2.x) && i1.x <= max(pa1.x, pa2.x) &&
                 i1.y >= min(pa1.y, pa2.y) && i1.y <= max(pa1.y, pa2.y)) &&
                (i2.x >= min(pa1.x, pa2.x) && i2.x <= max(pa1.x, pa2.x) &&
                 i2.y >= min(pa1.y, pa2.y) && i2.y <= max(pa1.y, pa2.y))) {
                return true;
            }
        }
    }

    return false;
}

MBR getMBR(vector<Point> puntos) {
    Point min, max;
    sort(puntos.begin(), puntos.end(), menorx);
    min.x = puntos.front().x;
    max.x = puntos.back().x;
    sort(puntos.begin(), puntos.end(), menory);
    min.y = puntos.front().y;
    max.y = puntos.back().y;
    return MBR(min, max);
}

MBR getMBRfromEdges(const vector<Edge>& edges) {
    vector<Point> puntos;
    for (const auto& edge : edges) {
        puntos.emplace_back(edge.start->lon, edge.start->lat);
        puntos.emplace_back(edge.end->lon, edge.end->lat);
    }
    return getMBR(puntos);
}

unordered_map<long long, Node> readNodes(const string& filename) {
    unordered_map<long long, Node> nodes;
    ifstream file(filename);
    string line;
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return nodes;
    }

    getline(file, line); // Skip header
    while (getline(file, line)) {
        istringstream ss(line);
        string token;
        Node node;

        getline(ss, token, ',');
        node.id = stoll(token);
        getline(ss, token, ',');
        node.lon = round(stold(token) * 1e7) / 1e7;
        getline(ss, token, ',');
        node.lat = round(stold(token) * 1e7) / 1e7;
        nodes[node.id] = node;
    }

    file.close();
    return nodes;
}

vector<Edge> readEdges(const string& filename, unordered_map<long long, Node>& nodes) {
    vector<Edge> edges;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return edges;
    }

    getline(file, line); 

    while (getline(file, line)) {
        istringstream ss(line);
        string token;
        Edge edge;

        try {
            getline(ss, token, ',');
            edge.id = stoll(token);
            getline(ss, token, ',');
            long long startId = stoll(token);
            getline(ss, token, ',');
            long long endId = stoll(token);

            edge.start = &nodes.at(startId);
            edge.end = &nodes.at(endId);
            edges.push_back(edge);
        }
        catch (const exception& e) {
            cout << "Error reading line: " << line << " - " << e.what() << endl;
        }
    }

    file.close();
    return edges;
}


struct Nodokd {
    Nodokd* izq, * der;
    MBR mbr;
    bool eshoja;
    Edge* inter;
    Nodokd(bool eshoja_ = false, Edge* inter_ = nullptr) : izq(nullptr), der(nullptr), eshoja(eshoja_), inter(inter_) {}
};

struct kd {
    Nodokd* raiz;

    void crearKD(vector<Edge>& aristas) {
        MBR mbr = getMBRfromEdges(aristas);
        raiz = construirKD(aristas, mbr, 0);
        cout << "kd construido" << endl;
    }

    Nodokd* construirKD(vector<Edge>& aristas, MBR mbr, int depth) {
        if (aristas.empty()) {
            return nullptr;
        }

        if (aristas.size() == 1) {
            Nodokd* nodo = new Nodokd(true, new Edge(aristas[0]));
            nodo->mbr = mbr;
            return nodo;
        }

        Nodokd* nodo = new Nodokd();
        nodo->mbr = mbr;

        vector<Edge> izqAristas, derAristas;
        MBR mbri, mbrd;

        if (depth % 2 == 0) {
            long double midX = (mbr.max.x + mbr.min.x) / 2;
            mbri = { mbr.min, {midX, mbr.max.y} };
            mbrd = { {midX, mbr.min.y}, mbr.max };
        }
        else {
            long double midY = (mbr.max.y + mbr.min.y) / 2;
            mbri = { {mbr.min.x, midY}, mbr.max };
            mbrd = { mbr.min, {mbr.max.x, midY} };
        }

        for (auto& arista : aristas) {
            if (Contenida(arista, mbri)) {
                izqAristas.push_back(arista);
            }
            if (Contenida(arista, mbrd)) {
                derAristas.push_back(arista);
            }
        }

        nodo->izq = construirKD(izqAristas, mbri, depth + 1);
        nodo->der = construirKD(derAristas, mbrd, depth + 1);

        return nodo;
    }

    Edge NN(Point q) {
        return AristaCercana(raiz, q, 0);
    }

    Edge AristaCercana(Nodokd* n, Point q, int depth, Nodokd* padre = nullptr, Nodokd* hermano = nullptr) {
        if (n == nullptr) {
            if (hermano != nullptr) {
                return AristaCercana(hermano, q, depth, padre, nullptr);
            } else if (padre != nullptr) {
                return AristaCercana(padre, q, depth - 1, nullptr, nullptr);
            } else {
                throw runtime_error("No se encontró una arista válida.");
            }
        }

        if (n->eshoja) {
            return *(n->inter);
        }

        Nodokd* next = nullptr;
        Nodokd* other = nullptr;

        if (depth % 2 == 0) {
            if (q.x < (n->mbr.max.x + n->mbr.min.x) / 2) {
                next = n->izq;
                other = n->der;
            } else {
                next = n->der;
                other = n->izq;
            }
        } else {
            if (q.y < (n->mbr.max.y + n->mbr.min.y) / 2) {
                next = n->izq;
                other = n->der;
            } else {
                next = n->der;
                other = n->izq;
            }
        }

        Edge best = AristaCercana(next, q, depth + 1, n, other);

        // Verificación si el MBR del otro nodo puede contener una mejor arista
        MBR mbrToCheck = (depth % 2 == 0) ?
            MBR({ (n->mbr.max.x + n->mbr.min.x) / 2, n->mbr.min.y }, { n->mbr.max.x, n->mbr.max.y }) :
            MBR({ n->mbr.min.x, (n->mbr.min.y + n->mbr.max.y) / 2 }, { n->mbr.max.x, n->mbr.max.y });

        if (puntoDentroDeMBR(q, mbrToCheck)) {
            try {
                Edge otherBest = AristaCercana(other, q, depth + 1, n, next);
                if (distancia(q, { otherBest.start->lon, otherBest.start->lat }) < distancia(q, { best.start->lon, best.start->lat })) {
                    best = otherBest;
                }
            } catch (const exception& e) {
                // Si no se encuentra una mejor arista en el hermano, continuamos con la mejor actual
            }
        }

        return best;
    }

    long double distancia(Point a, Point b) {
        return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    }

    void printKD() {
        if (!raiz) return;
        queue<Nodokd*> q;
        q.push(raiz);
        while (!q.empty()) {
            int levelSize = q.size();
            while (levelSize--) {
                Nodokd* nodo = q.front(); q.pop();
                cout << fixed << setprecision(7); 
                cout << "MBR: [(" << nodo->mbr.min.x << ", " << nodo->mbr.min.y << "), ("
                     << nodo->mbr.max.x << ", " << nodo->mbr.max.y << ")]";
                if (nodo->eshoja) {
                    cout << " -> Arista: [ " << nodo->inter->id << " (" << nodo->inter->start->lon << ", " << nodo->inter->start->lat << "), ("
                         << nodo->inter->end->lon << ", " << nodo->inter->end->lat << ")]";
                }
                cout << endl;
                if (nodo->izq) q.push(nodo->izq);
                if (nodo->der) q.push(nodo->der);
            }
            cout << "--------------------------------" << endl;
        }
    }
};

int main() {
    string nodesFile = "Recursos//nodes_reduced.csv";
    string edgesFile = "Recursos//edges_reduced.csv";

    unordered_map<long long, Node> nodes = readNodes(nodesFile); // Cambiado a long long
    vector<Edge> edges = readEdges(edgesFile, nodes);

    kd kdtree;
    kdtree.crearKD(edges);
    kdtree.printKD();

    Point q(-71.5631345, -16.4317946);
    Edge NN = kdtree.NN(q);

    cout << fixed << setprecision(7);
    cout << "query: (" << q.x << ", " << q.y << ") -> Arista mas cercana: " << NN.id << endl;

    return 0;
}

