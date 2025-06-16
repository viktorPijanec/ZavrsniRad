#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <getopt.h>
#include "bioparser/fasta_parser.hpp"
#include <fstream>
#include <iomanip>
#include <limits>

using namespace std;

namespace std
{
    template <>
    struct hash<std::pair<string, string>>
    {
        std::size_t operator()(const std::pair<string, string> &p) const noexcept
        {
            return std::hash<string>()(p.first) ^ (std::hash<string>()(p.second) << 1);
        }
    };

    template <>
    struct hash<std::pair<char, char>>
    {
        std::size_t operator()(const std::pair<char, char> &p) const noexcept
        {
            return std::hash<char>()(p.first) ^ (std::hash<char>()(p.second) << 1);
        }
    };

    template <>
    struct hash<std::pair<int, int>>
    {
        std::size_t operator()(const std::pair<int, int> &p) const noexcept
        {
            return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
        }
    };

    template <>
    struct hash<std::set<int>>
    {
        std::size_t operator()(const std::set<int> &s) const noexcept
        {
            std::size_t seed = 0;
            for (int i : s)
            {
                seed ^= std::hash<int>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

// constants
double GAP_OPEN = -10, GAP_EXT = -1;
bool NJ_tree = false;
string output_file = "izgraden_msa.txt", tree_file = "guide_tree.txt";

void print_help()
{
    cerr << "Usage: ./msa [options] <file>\n"
         << "Options:\n"
         << "  -h, --help       Show this help message and exit\n"
         << "  -g, --gap        Value for gap open penalty (default: -10)\n"
         << "  -e, --ext        Value for gap extend penalty (default: -1)\n"
         << "  -n, --njt        Use neighbor-joining tree instead of UPGMA\n"
         << "  -o, --out        Output file (default: izgraden_msa.txt)\n"
         << "  -t, --tree       Output tree (default: guide_tree.txt)\n"
         << "\n";
}

unordered_map<pair<char, char>, int> BLOSUM62;

void InitBlosum62()
{
    vector<string> residues = {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
    vector<vector<int>> scores = {
        // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},      // A
        {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},      // R
        {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},          // N
        {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},     // D
        {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},  // C
        {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},         // Q
        {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},        // E
        {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},    // G
        {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},      // H
        {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},     // I
        {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},     // L
        {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},      // K
        {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},      // M
        {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},      // F
        {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2}, // P
        {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2},         // S
        {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},     // T
        {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},  // W
        {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},    // Y
        {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4}       // V
    };

    for (size_t i = 0; i < residues.size(); ++i)
    {
        for (size_t j = 0; j < residues.size(); ++j)
        {
            char a = residues[i][0], b = residues[j][0];
            BLOSUM62[{a, b}] = scores[i][j];
            BLOSUM62[{b, a}] = scores[i][j]; // Ensure symmetry
        }
    }
}

// struct for bioparser
struct FastaSequence
{
private:
    const char *seq_name;
    const char *seq_data;
    std ::uint32_t name_len;
    std ::uint32_t data_len;

public:
    std::string name;
    std::string data;
    FastaSequence(const char *seq_name, std::uint32_t name_len, const char *seq_data, std::uint32_t data_len) : name(seq_name, name_len), data(seq_data, data_len), data_len(data_len) {}
    std::size_t length() const { return data_len; }
};

int CalcAlign(vector<char> &col1, vector<char> &col2)
{
    unordered_map<char, double> freq1, freq2;
    for (char c : col1)
        freq1[c] += 1.0;
    for (char c : col2)
        freq2[c] += 1.0;
    for (auto &[c1, f1] : freq1)
        f1 /= col1.size();
    for (auto &[c2, f2] : freq2)
        f2 /= col2.size();

    double score = 0;
    for (auto &[a, f1] : freq1)
    {
        for (auto &[b, f2] : freq2)
        {
            if (a == '-' || b == '-')
                score += f1 * f2 * GAP_EXT;
            else
                score += f1 * f2 * BLOSUM62[{a, b}];
        }
    }
    return static_cast<int>(score);
}

void AlignMultiple(vector<vector<char>> &seq1, vector<vector<char>> &seq2, vector<vector<char>> &aligned)
{
    enum DIRECTION
    {
        UP,
        LEFT,
        DIJAGONAL
    };
    struct Cell
    {
        int score;
        DIRECTION direction;
    };

    vector<vector<Cell>> matrica(seq1.size() + 1, vector<Cell>(seq2.size() + 1));

    matrica[0][0].score = GAP_OPEN;
    for (int i = 1; i <= seq2.size(); ++i)
    {
        matrica[0][i].score = matrica[0][i - 1].score + GAP_EXT;
        matrica[0][i].direction = LEFT;
    }
    for (int i = 1; i <= seq1.size(); ++i)
    {
        matrica[i][0].score = matrica[i - 1][0].score + GAP_EXT;
        matrica[i][0].direction = UP;
    }

    for (int i = 1; i <= seq1.size(); ++i)
    {
        for (int j = 1; j <= seq2.size(); ++j)
        {
            int gore, lijevo;
            if (matrica[i - 1][j].direction == UP || matrica[i - 1][j].direction == LEFT)
                gore = matrica[i - 1][j].score + GAP_EXT;
            else
                gore = matrica[i - 1][j].score + GAP_OPEN;
            if (matrica[i][j - 1].direction == UP || matrica[i][j - 1].direction == LEFT)
                lijevo = matrica[i][j - 1].score + GAP_EXT;
            else
                lijevo = matrica[i][j - 1].score + GAP_OPEN;
            int dijagonala = matrica[i - 1][j - 1].score + CalcAlign(seq1[i - 1], seq2[j - 1]);
            if (max(max(gore, lijevo), dijagonala) == gore)
            {
                matrica[i][j].score = gore;
                matrica[i][j].direction = UP;
            }
            else if (max(max(gore, lijevo), dijagonala) == lijevo)
            {
                matrica[i][j].score = lijevo;
                matrica[i][j].direction = LEFT;
            }
            else
            {
                matrica[i][j].score = dijagonala;
                matrica[i][j].direction = DIJAGONAL;
            }
        }
    }
    // if (seq1[0].size() == 1 && seq2[0].size() == 1)
    // {
    //     cout << "matrica: " << endl;
    //     for (int i = 0; i <= seq1.size(); ++i)
    //     {
    //         for (int j = 0; j <= seq2.size(); ++j)
    //         {
    //             cout << matrica[i][j].score << " ";
    //         }
    //         cout << endl;
    //     }
    // }

    // finding best alignment - used for semi global
    // int maxalign = 0, maxi = 0, maxj = 0;
    // for (int i = 0; i <= seq1.size(); ++i)
    // {
    //     if (matrica[i][seq2.size()].score > maxalign)
    //     {
    //         maxalign = matrica[i][seq2.size()].score;
    //         maxi = i;
    //         maxj = seq2.size();
    //     }
    // }
    // for (int i = 0; i <= seq2.size(); ++i)
    // {
    //     if (matrica[seq1.size()][i].score > maxalign)
    //     {
    //         maxalign = matrica[seq1.size()][i].score;
    //         maxi = seq1.size();
    //         maxj = i;
    //     }
    // }

    // building alignment
    int treni = seq1.size(), trenj = seq2.size();
    while (treni > 0 || trenj > 0)
    {
        if (matrica[treni][trenj].direction == DIJAGONAL)
        {
            vector<char> temp = seq1[treni - 1];
            temp.insert(temp.end(), seq2[trenj - 1].begin(), seq2[trenj - 1].end());
            aligned.insert(aligned.begin(), temp);
            --treni;
            --trenj;
        }
        else if (matrica[treni][trenj].direction == UP)
        {
            vector<char> temp = seq1[treni - 1];
            for (int i = 0; i < seq2[0].size(); ++i)
                temp.push_back('-');
            aligned.insert(aligned.begin(), temp);
            --treni;
        }
        else if (matrica[treni][trenj].direction == LEFT)
        {
            vector<char> temp;
            for (int i = 0; i < seq1[0].size(); ++i)
                temp.push_back('-');
            temp.insert(temp.end(), seq2[trenj - 1].begin(), seq2[trenj - 1].end());
            aligned.insert(aligned.begin(), temp);
            --trenj;
        }
    }
}

int Align(string &seq1, string &seq2, string *path)
{
    enum DIRECTION
    {
        UP,
        LEFT,
        DIJAGONAL
    };
    struct Cell
    {
        int score;
        DIRECTION direction;
    };

    vector<vector<Cell>> matrica(seq1.size() + 1, vector<Cell>(seq2.size() + 1));

    matrica[0][0].score = GAP_OPEN;
    for (int i = 1; i <= seq2.size(); ++i)
    {
        matrica[0][i].score = matrica[0][i - 1].score + GAP_EXT;
        matrica[0][i].direction = LEFT;
    }
    for (int i = 1; i <= seq1.size(); ++i)
    {
        matrica[i][0].score = matrica[i - 1][0].score + GAP_EXT;
        matrica[i][0].direction = UP;
    }

    for (int i = 1; i <= seq1.size(); ++i)
    {
        for (int j = 1; j <= seq2.size(); ++j)
        {
            int gore, lijevo;
            if (matrica[i - 1][j].direction == UP || matrica[i - 1][j].direction == LEFT)
                gore = matrica[i - 1][j].score + GAP_EXT;
            else
                gore = matrica[i - 1][j].score + GAP_OPEN;
            if (matrica[i][j - 1].direction == UP || matrica[i][j - 1].direction == LEFT)
                lijevo = matrica[i][j - 1].score + GAP_EXT;
            else
                lijevo = matrica[i][j - 1].score + GAP_OPEN;
            char a = seq1[i - 1];
            char b = seq2[j - 1];
            int sub_score = (a == '-' || b == '-') ? GAP_EXT : BLOSUM62[{a, b}];
            int dijagonala = matrica[i - 1][j - 1].score + sub_score;
            if (max(max(gore, lijevo), dijagonala) == gore)
            {
                matrica[i][j].score = gore;
                matrica[i][j].direction = UP;
            }
            else if (max(max(gore, lijevo), dijagonala) == lijevo)
            {
                matrica[i][j].score = lijevo;
                matrica[i][j].direction = LEFT;
            }
            else
            {
                matrica[i][j].score = dijagonala;
                matrica[i][j].direction = DIJAGONAL;
            }
        }
    }
    /*
    if (seq1[0] == 'M' && seq2[0] == 'M' && seq1[1] == 'Q' && seq2[1] == 'H')
    {
        cout << "matrica: " << endl;
        for (int i = 0; i <= seq1.size(); ++i)
        {
            for (int j = 0; j <= seq2.size(); ++j)
            {
                cout << matrica[i][j].score << " ";
            }
            cout << endl;
        }
    }*/

    if (path != nullptr)
    {
        int treni = seq1.size(), trenj = seq2.size();
        while (treni > 0 || trenj > 0)
        {
            if (matrica[treni][trenj].direction == DIJAGONAL)
            {
                path->insert(0, "D");
                treni--;
                trenj--;
            }
            else if (matrica[treni][trenj].direction == UP)
            {
                path->insert(0, "U");
                treni--;
            }
            else
            {
                path->insert(0, "L");
                trenj--;
            }
        }
    }

    return matrica[seq1.size()][seq2.size()].score;
}

bool PairContainsString(pair<string, string> p, string str)
{
    if (p.first == str || p.second == str)
        return true;
    return false;
}

string BuildTree(unordered_map<pair<string, string>, double> &initial_alignments)
{
    unordered_map<pair<string, string>, double> alignments = initial_alignments;
    set<string> set_svih_cvorova;

    for (const auto &[key, value] : alignments)
    {
        set_svih_cvorova.insert(key.first);
        set_svih_cvorova.insert(key.second);
    }

    unordered_map<string, set<int>> atomi;
    unordered_map<string, string> ime_do_stringa;

    string ret;
    while (!alignments.empty())
    {
        // trazimo max
        double max_val = alignments.begin()->second;
        pair<string, string> max_pair = alignments.begin()->first;
        for (const auto &[key, value] : alignments)
        {
            if (value > max_val)
            {
                max_val = value;
                max_pair = key;
            }
        }

        // cout << "max pair: " << max_pair.first << " " << max_pair.second << endl;

        // brisanje para iz mape
        alignments.erase(max_pair);

        // brisanje cvora iz seta svih cvorova
        set_svih_cvorova.erase(max_pair.first);
        set_svih_cvorova.erase(max_pair.second);

        string ime_novog_cvora = "";
        ime_novog_cvora += max_pair.first + "_" + max_pair.second;
        set<int> clanovi_novog_cvora;
        string str_novi_cvor = "(";
        // provjera jel prvi clan atom ili cvor
        if (max_pair.first.find('_') == string::npos)
        {
            clanovi_novog_cvora.insert(stoi(max_pair.first));
            str_novi_cvor += max_pair.first;
        }
        else
        {
            clanovi_novog_cvora.insert(atomi[max_pair.first].begin(), atomi[max_pair.first].end());
            str_novi_cvor += ime_do_stringa[max_pair.first];
        }
        str_novi_cvor += ",";
        // provjera jel drugi clan atom ili cvor
        if (max_pair.second.find('_') == string::npos)
        {
            clanovi_novog_cvora.insert(stoi(max_pair.second));
            str_novi_cvor += max_pair.second;
        }
        else
        {
            clanovi_novog_cvora.insert(atomi[max_pair.second].begin(), atomi[max_pair.second].end());
            str_novi_cvor += ime_do_stringa[max_pair.second];
        }
        str_novi_cvor += ")\n";
        // dodavanje clanova u mapu
        atomi[ime_novog_cvora] = clanovi_novog_cvora;
        ime_do_stringa[ime_novog_cvora] = str_novi_cvor;
        ret = str_novi_cvor;
        // brisanje i stvaranje novih parova kao srednjih vrijednosti
        for (auto str : set_svih_cvorova)
        {
            double avg = 0;
            int brojac = 0;
            // cout << "tren str: " << str << endl;
            if (str.find('_') == string::npos)
            {
                if (max_pair.first.find('_') == string::npos)
                {
                    avg += initial_alignments[make_pair(min(str, max_pair.first), max(str, max_pair.first))];
                    brojac++;
                }
                else
                {
                    for (int i : atomi[max_pair.first])
                    {
                        avg += initial_alignments[make_pair(min(str, to_string(i)), max(str, to_string(i)))];
                        brojac++;
                    }
                }
                if (max_pair.second.find('_') == string::npos)
                {
                    avg += initial_alignments[make_pair(min(str, max_pair.second), max(str, max_pair.second))];
                    brojac++;
                }
                else
                {
                    for (int i : atomi[max_pair.second])
                    {
                        avg += initial_alignments[make_pair(min(str, to_string(i)), max(str, to_string(i)))];
                        brojac++;
                    }
                }
            }
            else
            {
                for (int i : atomi[str])
                {
                    if (max_pair.first.find('_') == string::npos)
                    {
                        avg += initial_alignments[make_pair(min(to_string(i), max_pair.first), max(to_string(i), max_pair.first))];
                        brojac++;
                    }
                    else
                    {
                        for (int j : atomi[max_pair.first])
                        {
                            avg += initial_alignments[make_pair(to_string(min(i, j)), to_string(max(i, j)))];
                            brojac++;
                        }
                    }
                    if (max_pair.second.find('_') == string::npos)
                    {
                        avg += initial_alignments[make_pair(min(to_string(i), max_pair.second), max(to_string(i), max_pair.second))];
                        brojac++;
                    }
                    else
                    {
                        for (int j : atomi[max_pair.second])
                        {
                            avg += initial_alignments[make_pair(to_string(min(i, j)), to_string(max(i, j)))];
                            brojac++;
                        }
                    }
                }
            }
            if (str.find('_') == string::npos && max_pair.first.find('_') == string::npos)
                alignments.erase(make_pair(to_string(min(stoi(str), stoi(max_pair.first))), to_string(max(stoi(str), stoi(max_pair.first)))));
            else
                alignments.erase(make_pair(min(str, max_pair.first), max(str, max_pair.first)));
            if (str.find('_') == string::npos && max_pair.second.find('_') == string::npos)
                alignments.erase(make_pair(to_string(min(stoi(str), stoi(max_pair.second))), to_string(max(stoi(str), stoi(max_pair.second)))));
            else
                alignments.erase(make_pair(min(str, max_pair.second), max(str, max_pair.second)));
            avg /= brojac;
            alignments[make_pair(min(ime_novog_cvora, str), max(ime_novog_cvora, str))] = avg;
        }
        // for (const auto &[key, value] : alignments)
        // {
        //     cout << key.first << " " << key.second << " " << value << endl;
        // }
        set_svih_cvorova.insert(ime_novog_cvora);
    }
    return ret;
}

string NeighborJoiningTree(const unordered_map<pair<int, int>, double> &distances, int num_seqs)
{
    vector<set<int>> clusters;
    unordered_map<int, string> names;
    unordered_map<pair<int, int>, double> D = distances;

    for (int i = 0; i < num_seqs; ++i)
    {
        clusters.push_back({i});
        names[i] = to_string(i);
    }

    int next_id = num_seqs;

    while (clusters.size() > 2)
    {
        int n = clusters.size();
        vector<int> ids;
        for (const auto &cl : clusters)
            ids.push_back(*cl.begin());

        // Compute Q-matrix
        unordered_map<pair<int, int>, double> Q;
        unordered_map<int, double> total;
        for (int i : ids)
            for (int j : ids)
                if (i != j)
                    total[i] += D[{min(i, j), max(i, j)}];

        double min_Q = numeric_limits<double>::infinity();
        int i_min = -1, j_min = -1;

        for (int i = 0; i < ids.size(); ++i)
        {
            for (int j = i + 1; j < ids.size(); ++j)
            {
                int a = ids[i], b = ids[j];
                double q = (n - 2) * D[{min(a, b), max(a, b)}] - total[a] - total[b];
                Q[{a, b}] = q;
                if (q < min_Q)
                {
                    min_Q = q;
                    i_min = a;
                    j_min = b;
                }
            }
        }

        // Merge clusters i_min and j_min
        string new_name = "(" + names[i_min] + "," + names[j_min] + ")";
        names[next_id] = new_name;

        // Update distances
        for (int k : ids)
        {
            if (k == i_min || k == j_min)
                continue;
            double d1 = D[{min(i_min, k), max(i_min, k)}];
            double d2 = D[{min(j_min, k), max(j_min, k)}];
            double d_new = 0.5 * (d1 + d2 - D[{min(i_min, j_min), max(i_min, j_min)}]);
            D[{min(next_id, k), max(next_id, k)}] = d_new;
        }

        // Erase old distances
        for (int k : ids)
        {
            D.erase({min(i_min, k), max(i_min, k)});
            D.erase({min(j_min, k), max(j_min, k)});
        }
        D.erase({min(i_min, j_min), max(i_min, j_min)});

        // Update cluster list
        set<int> new_cluster = {next_id};
        clusters.erase(remove_if(clusters.begin(), clusters.end(), [&](const set<int> &c)
                                 { return c.count(i_min) || c.count(j_min); }),
                       clusters.end());
        clusters.push_back(new_cluster);

        ++next_id;
    }

    // Join the last two clusters
    int a = *clusters[0].begin(), b = *clusters[1].begin();
    return "(" + names[a] + "," + names[b] + ")";
}

int CountInString(string str, char c)
{
    int ret = 0;
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == c)
            ret++;
    }
    return ret;
}

vector<string> SplitString(string str, char c)
{
    vector<string> ret;
    string temp = "";
    for (int i = 0; i < str.size(); i++)
    {
        if (str[i] == c)
        {
            ret.push_back(temp);
            temp = "";
        }
        else
            temp += str[i];
    }
    ret.push_back(temp);
    return ret;
}

string BuildByTree(string stablo, unordered_map<string, vector<vector<char>>> &filogen_stablo_map, string &zadnji, const vector<unique_ptr<FastaSequence>> &sequences)
{
    if (CountInString(stablo, ',') == 1)
    {
        vector<string> podijela = SplitString(stablo, ',');
        string prvi = podijela[0].substr(1), drugi = podijela[1].substr(0, podijela[1].size() - 1);
        zadnji = prvi + '_' + drugi;
        if (filogen_stablo_map.find(prvi) == filogen_stablo_map.end())
        {
            vector<vector<char>> temp1;
            for (int i = 0; i < sequences[stoi(prvi)]->data.size(); i++)
            {
                vector<char> temp2;
                temp2.push_back(sequences[stoi(prvi)]->data[i]);
                temp1.push_back(temp2);
            }
            filogen_stablo_map[prvi] = temp1;
        }
        if (filogen_stablo_map.find(drugi) == filogen_stablo_map.end())
        {
            vector<vector<char>> temp1;
            for (int i = 0; i < sequences[stoi(drugi)]->data.size(); i++)
            {
                vector<char> temp2;
                temp2.push_back(sequences[stoi(drugi)]->data[i]);
                temp1.push_back(temp2);
            }
            filogen_stablo_map[drugi] = temp1;
        }
        vector<vector<char>> temp;
        AlignMultiple(filogen_stablo_map[prvi], filogen_stablo_map[drugi], temp);
        // ispis tempa
        // for (int i = 0; i < temp[0].size(); i++)
        // {
        //     for (int j = 0; j < temp.size(); j++)
        //     {
        //         cout << temp[j][i];
        //     }
        //     cout << endl;
        // }
        filogen_stablo_map[zadnji] = temp;
        return zadnji;
    }
    else
    {
        int l_zagrada = 0;
        int d_zagrada = 0;
        for (int i = 0; i < stablo.size(); i++)
        {
            if (stablo[i] == ',' && l_zagrada - d_zagrada == 1)
            {
                string temp1 = BuildByTree(stablo.substr(1, i - 1), filogen_stablo_map, zadnji, sequences);
                string temp2 = BuildByTree(stablo.substr(i + 1, stablo.size() - i - 2), filogen_stablo_map, zadnji, sequences);
                stablo = '(' + temp1 + ',' + temp2 + ')';
                stablo = BuildByTree(stablo, filogen_stablo_map, zadnji, sequences);
            }
            else if (stablo[i] == '(')
                l_zagrada++;
            else if (stablo[i] == ')')
                d_zagrada++;
        }
        return stablo;
    }
}

int main(int argc, char *argv[])
{
    // Define long options
    const struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"gap", required_argument, nullptr, 'g'},
        {"ext", required_argument, nullptr, 'e'},
        {"njt", no_argument, nullptr, 'n'},
        {"out", required_argument, nullptr, 'o'},
        {"tree", required_argument, nullptr, 't'},
        {nullptr, 0, nullptr, 0}};

    int option_index = 0;
    int opt;

    // Parse command line arguments
    while ((opt = getopt_long(argc, argv, "hg:e:o:t:n", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'h':
            print_help();
            return 0;
        case 'g':
            GAP_OPEN = atoi(optarg);
            break;
        case 'e':
            GAP_EXT = atoi(optarg);
            break;
        case 'n':
            NJ_tree = true;
            break;
        case 'o':
            output_file = optarg;
            break;
        case 't':
            tree_file = optarg;
            break;
        case '?':
            print_help();
            return 1;
        default:
            abort();
        }
    }

    // After options, we expect two additional arguments (reference_file and fragments_file)
    if (argc - optind < 1)
    {
        cerr << "Error: Missing required file arguments.\n";
        print_help();
        return 1;
    }

    std::ios_base::sync_with_stdio(false);

    InitBlosum62();

    // creating parser
    auto parser = bioparser::Parser<FastaSequence>::Create<bioparser::FastaParser>(argv[argc - 1]);

    // parse file
    auto sequences = parser->Parse(-1);

    /*
    for (int i = 0; i < sequences.size(); i++)
    {
        cout << sequences[i]->name << endl
             << sequences[i]->data << endl
             << endl;
    }
    */

    unordered_map<pair<string, string>, double> alignments;
    double max_align = numeric_limits<double>::min();
    for (int i = 0; i < sequences.size(); ++i)
    {
        for (int j = i + 1; j < sequences.size(); ++j)
        {
            pair<string, string> par = make_pair(to_string(i), to_string(j));
            string trenpath;
            alignments.insert({par, (double)Align(sequences[i]->data, sequences[j]->data, &trenpath)});
            if ((double)alignments[par] > max_align)
                max_align = (double)alignments[par];
            // cout << "Align(" << i + 1 << "," << j + 1 << ") : " << alignments[par] << endl;
        }
    }

    cout << "Building filogenic tree..." << endl;
    string filogen_stablo;
    if (NJ_tree)
    {
        // creating distance matrix from similarity
        unordered_map<pair<int, int>, double> distance_matrix;
        for (auto &[pair, score] : alignments)
        {
            int i = stoi(pair.first), j = stoi(pair.second);
            distance_matrix[{min(i, j), max(i, j)}] = max_align - score;
        }

        filogen_stablo = NeighborJoiningTree(distance_matrix, sequences.size());
    }
    else
    {
        filogen_stablo = BuildTree(alignments);
    }
    vector<int> final_order;
    int poc_ind_broja = -1, broj_pronaden = 0;
    for (int i = 0; i < filogen_stablo.size(); i++)
    {
        if (filogen_stablo[i] == ',' || filogen_stablo[i] == ')' || filogen_stablo[i] == '(' || filogen_stablo[i] == '\n')
        {
            broj_pronaden = 0;
            continue;
        }
        if (broj_pronaden == 0)
        {
            poc_ind_broja = i;
            broj_pronaden = 1;
        }
        if (i != filogen_stablo.size() - 1 && filogen_stablo[i + 1] >= '0' && filogen_stablo[i + 1] <= '9')
            continue;
        final_order.push_back(stoi(filogen_stablo.substr(poc_ind_broja, i - poc_ind_broja + 1)));
    }
    /*
    for (int i = 0; i < final_order.size(); i++)
    {
        cout << final_order[i] << " ";
    }
    cout << endl;
    */
    ofstream filogen(tree_file);
    filogen << filogen_stablo;
    filogen.close();

    unordered_map<string, vector<vector<char>>> filogen_stablo_map;
    string zadnji_kljuc;

    cout << "Building by tree..." << endl;
    BuildByTree(filogen_stablo, filogen_stablo_map, zadnji_kljuc, sequences);

    ofstream complete_msa(output_file);

    complete_msa << "PileUp" << endl
                 << endl
                 << "   MSF:" << std::setw(5) << std::right << filogen_stablo_map[zadnji_kljuc].size()
                 << "  Type: P    Check: " << std::setw(5) << std::right << 0 << "   .." << endl
                 << endl;

    for (int i = 0; i < sequences.size(); i++)
    {
        complete_msa << " " << "Name: " << sequences[i]->name << " oo  Len: " << filogen_stablo_map[zadnji_kljuc].size() << "  Check: " << std::setw(5) << std::right << 0 << "  Weight:  10.0" << endl;
    }
    complete_msa << endl
                 << "//" << endl
                 << endl;

    int izadi_van = 0;
    int okreti = 0;
    while (!izadi_van)
    {
        for (int i = 0; i < filogen_stablo_map[zadnji_kljuc][0].size(); i++)
        {
            complete_msa << std::setw(12) << std::left << sequences[final_order[i]]->name;
            for (int j = 0; j < 50; j++)
            {
                if (j % 10 == 0 && j != 0)
                    complete_msa << " ";
                if (okreti * 50 + j >= filogen_stablo_map[zadnji_kljuc].size())
                {
                    izadi_van = 1;
                    break;
                }
                complete_msa << filogen_stablo_map[zadnji_kljuc][okreti * 50 + j][i];
            }
            complete_msa << endl;
        }
        complete_msa << endl
                     << endl;
        okreti++;
    }
    complete_msa.close();

    cout << "output saved in " << output_file << endl;

    return 0;
}