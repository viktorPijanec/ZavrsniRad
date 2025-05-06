#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "bioparser/fasta_parser.hpp"

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
int GAP = -1, MATCH = 1, MISMATCH = 0;

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

int CalcAlign(vector<char> &seq1, vector<char> &seq2)
{
    int ret = 0;
    for (int i = 0; i < seq1.size(); i++)
    {
        for (int j = i + 1; j < seq1.size(); j++)
        {
            if (seq1[i] == seq1[j])
                ret += MATCH;
            else
                ret += MISMATCH;
        }
        for (int j = 0; j < seq2.size(); j++)
        {
            if (seq1[i] == seq2[j])
                ret += MATCH;
            else
                ret += MISMATCH;
        }
    }
    for (int i = 0; i < seq2.size(); i++)
    {
        for (int j = i + 1; j < seq2.size(); j++)
        {
            if (seq2[i] == seq2[j])
                ret += MATCH;
            else
                ret += MISMATCH;
        }
    }
    return ret;
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

    Cell matrica[seq1.size() + 1][seq2.size() + 1];

    for (int i = 0; i <= seq2.size(); ++i)
    {
        matrica[0][i].score = i * GAP;
        matrica[0][i].direction = LEFT;
    }
    for (int i = 1; i <= seq1.size(); ++i)
    {
        matrica[i][0].score = i * GAP;
        matrica[i][0].direction = UP;
    }

    for (int i = 1; i <= seq1.size(); ++i)
    {
        for (int j = 1; j <= seq2.size(); ++j)
        {
            int gore = matrica[i - 1][j].score + GAP, lijevo = matrica[i][j - 1].score + GAP;
            int dijagonala = CalcAlign(seq1[i - 1], seq2[j - 1]);
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

    Cell matrica[seq1.size() + 1][seq2.size() + 1];

    for (int i = 0; i <= seq2.size(); ++i)
    {
        matrica[0][i].score = i * GAP;
        matrica[0][i].direction = LEFT;
    }
    for (int i = 1; i <= seq1.size(); ++i)
    {
        matrica[i][0].score = i * GAP;
        matrica[i][0].direction = UP;
    }

    for (int i = 1; i <= seq1.size(); ++i)
    {
        for (int j = 1; j <= seq2.size(); ++j)
        {
            int gore = matrica[i - 1][j].score + GAP, lijevo = matrica[i][j - 1].score + GAP;
            int dijagonala = seq1[i - 1] == seq2[j - 1] ? matrica[i - 1][j - 1].score + MATCH : matrica[i - 1][j - 1].score + MISMATCH;
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

string BuildTree(unordered_map<pair<string, string>, double> initial_alignments)
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
        str_novi_cvor += ")";
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
            alignments.erase(make_pair(min(str, max_pair.first), max(str, max_pair.first)));
            alignments.erase(make_pair(min(str, max_pair.second), max(str, max_pair.second)));
            avg /= brojac;
            alignments[make_pair(min(ime_novog_cvora, str), max(ime_novog_cvora, str))] = avg;
            // for (const auto& [key,value] : alignments){
            //     cout << key.first << " " << key.second << " " << value << endl;
            // }
        }
        set_svih_cvorova.insert(ime_novog_cvora);
    }
    return ret;
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

string IzgradiPoStablu(string stablo, unordered_map<string, vector<vector<char>>> &filogen_stablo_map, string &zadnji, const vector<unique_ptr<FastaSequence>> &sequences)
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
        for (int i = 0; i < temp[0].size(); i++)
        {
            for (int j = 0; j < temp.size(); j++)
            {
                cout << temp[j][i];
            }
            cout << endl;
        }
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
                string temp1 = IzgradiPoStablu(stablo.substr(1, i - 1), filogen_stablo_map, zadnji, sequences);
                string temp2 = IzgradiPoStablu(stablo.substr(i + 1, stablo.size() - i - 2), filogen_stablo_map, zadnji, sequences);
                stablo = '(' + temp1 + ',' + temp2 + ')';
                stablo = IzgradiPoStablu(stablo, filogen_stablo_map, zadnji, sequences);
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
    // creating parser
    auto parser = bioparser::Parser<FastaSequence>::Create<bioparser::FastaParser>(argv[1]);

    // parse file
    auto sequences = parser->Parse(-1);

    for (int i = 0; i < sequences.size(); i++)
    {
        cout << sequences[i]->name << endl
             << sequences[i]->data << endl
             << endl;
    }

    unordered_map<pair<string, string>, double> alignments;

    int maxi, maxj, maxalign;
    string maxpath;
    for (int i = 0; i < sequences.size(); ++i)
    {
        for (int j = i + 1; j < sequences.size(); ++j)
        {
            pair<string, string> par = make_pair(to_string(i), to_string(j));
            string trenpath;
            alignments.insert({par, (double)Align(sequences[i]->data, sequences[j]->data, &trenpath)});
            if (i == 0 && j == 1)
            {
                maxi = i;
                maxj = j;
                maxalign = alignments[par];
                maxpath = trenpath;
            }
            else if (alignments[par] > maxalign)
            {
                maxi = i;
                maxj = j;
                maxalign = alignments[par];
                maxpath = trenpath;
            }
            cout << "Align(" << i + 1 << "," << j + 1 << ") : " << alignments[par] << endl;
        }
    }

    cout << "Building filogenic tree..." << endl;
    string filogen_stablo = BuildTree(alignments);
    cout << "Tree created." << endl;

    cout << filogen_stablo << endl;

    unordered_map<string, vector<vector<char>>> filogen_stablo_map;
    string zadnji_kljuc;

    IzgradiPoStablu(filogen_stablo, filogen_stablo_map, zadnji_kljuc, sequences);

    cout << "Izgraden msa:  " << endl
         << endl;
    for (int i = 0; i < filogen_stablo_map[zadnji_kljuc][0].size(); i++)
    {
        for (int j = 0; j < filogen_stablo_map[zadnji_kljuc].size(); j++)
        {
            cout << filogen_stablo_map[zadnji_kljuc][j][i];
        }
        cout << endl;
    }

    cout << endl
         << "Pokazni primjerak na najboljem paru (u ovom slucaju " + to_string(maxi + 1) + " i " + to_string(maxj + 1) + "):" << endl
         << endl;
    string novi_str1 = "", novi_str2 = "";
    int ind1 = 0, ind2 = 0;
    for (int i = 0; i < maxpath.length(); ++i)
    {
        if (maxpath[i] == 'D')
        {
            novi_str1 += sequences[maxi]->data[ind1++];
            novi_str2 += sequences[maxj]->data[ind2++];
        }
        else if (maxpath[i] == 'U')
        {
            novi_str1 += sequences[maxi]->data[ind1++];
            novi_str2 += '-';
        }
        else
        {
            novi_str1 += '-';
            novi_str2 += sequences[maxj]->data[ind2++];
        }
    }
    cout << novi_str1 << endl;
    cout << novi_str2 << endl;
    return 0;
}