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
            for (int i : s){
                seed ^= std::hash<int>()(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }
            return seed;
        }
    };
}

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

int Align(string seq1, string seq2, string *path)
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

    Cell matrica[seq1.length() + 1][seq2.length() + 1];

    for (int i = 0; i <= seq2.length(); ++i)
    {
        matrica[0][i].score = i * GAP;
        matrica[0][i].direction = LEFT;
    }
    for (int i = 1; i <= seq1.length(); ++i)
    {
        matrica[i][0].score = i * GAP;
        matrica[i][0].direction = UP;
    }

    for (int i = 1; i <= seq1.length(); ++i)
    {
        for (int j = 1; j <= seq2.length(); ++j)
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
        // cout << endl;
        // for (int i = 0; i < 10; ++i)
        // {
        //     for (int j = 0; j < 10; ++j)
        //     {
        //         printf("%3d ", matrica[i][j].score);
        //     }
        //     cout << '\n';
        // }
        int treni = seq1.length(), trenj = seq2.length();
        while (treni > 0 && trenj > 0)
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

    return matrica[seq1.length()][seq2.length()].score;
}

bool PairContainsString(pair<string,string> p, string str){
    if (p.first == str || p.second == str) return true;
    return false;
}

// string smallString(string str1, string str2){
//     if (str1<str2) return str1;
//     else return str2;
// }

// string bigString(string str1, string str2){
//     if (str1>str2) return str1;
//     else return str2;
// }

string BuildTree(unordered_map<pair<string,string>,double> initial_alignments)
{
    unordered_map<pair<string,string>,double> alignments = initial_alignments;
    set<string> set_svih_cvorova;

    for (const auto &[key,value] : alignments){
        set_svih_cvorova.insert(key.first);
        set_svih_cvorova.insert(key.second);
    }

    unordered_map<string, set<int>> atomi;
    unordered_map<string, string> ime_do_stringa;

    string ret;
    while(!alignments.empty()){
        // trazimo max
        double max_val = alignments.begin()->second;
        pair<string,string> max_pair = alignments.begin()->first;
        for (const auto &[key,value] : alignments){
            if (value > max_val){
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
        if (max_pair.first.find('_') == string::npos){
            clanovi_novog_cvora.insert(stoi(max_pair.first));
            str_novi_cvor += max_pair.first;
        }
        else{
            clanovi_novog_cvora.insert(atomi[max_pair.first].begin(),atomi[max_pair.first].end()); 
            str_novi_cvor += ime_do_stringa[max_pair.first];
        }
        str_novi_cvor += ",";
        if (max_pair.second.find('_') == string::npos){
            clanovi_novog_cvora.insert(stoi(max_pair.second));
            str_novi_cvor += max_pair.second;
        }
        else{
            clanovi_novog_cvora.insert(atomi[max_pair.second].begin(),atomi[max_pair.second].end()); 
            str_novi_cvor += ime_do_stringa[max_pair.second];
        }
        str_novi_cvor += ")";
        // dodavanje clanova u mapu
        atomi[ime_novog_cvora] = clanovi_novog_cvora;
        ime_do_stringa[ime_novog_cvora] = str_novi_cvor;
        ret = str_novi_cvor;
        // brisanje i stvaranje novih parova kao srednjih vrijednosti
        for (auto str : set_svih_cvorova){
            double avg = 0;
            int brojac = 0;
            // cout << "tren str: " << str << endl;
            if (str.find('_') == string::npos){
                if (max_pair.first.find('_') == string::npos){
                    avg += initial_alignments[make_pair(min(str,max_pair.first),max(str,max_pair.first))];
                    brojac++;
                }
                else{
                    for (int i : atomi[max_pair.first]){
                        avg+=initial_alignments[make_pair(min(str,to_string(i)),max(str,to_string(i)))];
                        brojac++;
                    }
                }
                if (max_pair.second.find('_') == string::npos){
                    avg += initial_alignments[make_pair(min(str,max_pair.second),max(str,max_pair.second))];
                    brojac++;
                }
                else{
                    for (int i : atomi[max_pair.second]){
                        avg+=initial_alignments[make_pair(min(str,to_string(i)),max(str,to_string(i)))];
                        brojac++;
                    }
                }
            }
            else{
                for (int i : atomi[str]){
                    if (max_pair.first.find('_') == string::npos){
                        avg += initial_alignments[make_pair(min(to_string(i),max_pair.first),max(to_string(i),max_pair.first))];
                        brojac++;
                    }
                    else{
                        for (int j : atomi[max_pair.first]){
                            avg += initial_alignments[make_pair(to_string(min(i,j)),to_string(max(i,j)))];
                            brojac++;
                        }
                    }
                    if (max_pair.second.find('_') == string::npos){
                        avg += initial_alignments[make_pair(min(to_string(i),max_pair.second),max(to_string(i),max_pair.second))];
                        brojac++;
                    }
                    else{
                        for (int j : atomi[max_pair.second]){
                            avg += initial_alignments[make_pair(to_string(min(i,j)),to_string(max(i,j)))];
                            brojac++;
                        }
                    }
                }
            }
            alignments.erase(make_pair(min(str,max_pair.first),max(str,max_pair.first)));
            alignments.erase(make_pair(min(str,max_pair.second),max(str,max_pair.second)));
            avg/=brojac;
            alignments[make_pair(min(ime_novog_cvora,str),max(ime_novog_cvora,str))] = avg;
            // for (const auto& [key,value] : alignments){
            //     cout << key.first << " " << key.second << " " << value << endl;
            // }
        }
        set_svih_cvorova.insert(ime_novog_cvora);
    }
    return ret;
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
    for (int i = 0; i < sequences.size(); ++i)
    {
        for (int j = i + 1; j < sequences.size(); ++j)
        {
            pair<string, string> par = make_pair(to_string(i), to_string(j));
            alignments.insert({par, (double)Align(sequences[i]->data, sequences[j]->data, nullptr)});
            if (i == 0 && j == 1)
            {
                maxi = i;
                maxj = j;
                maxalign = alignments[par];
            }
            else if (alignments[par] > maxalign)
            {
                maxi = i;
                maxj = j;
                maxalign = alignments[par];
            }
            cout << "Align(" << i + 1 << "," << j + 1 << ") : " << alignments[par] << endl;
        }
    }

    cout << "Building filogenic tree..." << endl;
    string filogen_stablo = BuildTree(alignments);
    cout << "Tree created." << endl;

    cout << filogen_stablo << endl;

    string pomstr;
    Align(sequences[maxi]->data, sequences[maxj]->data, &pomstr);
    cout << endl
         << "Pokazni primjerak na najboljem paru (u ovom slucaju " + to_string(maxi + 1) + " i " + to_string(maxj + 1) + "):" << endl
         << endl;
    // cout << pomstr << endl;
    string novi_str1 = "", novi_str2 = "";
    int ind1 = 0, ind2 = 0;
    for (int i = 0; i < pomstr.length(); ++i)
    {
        if (pomstr[i] == 'D')
        {
            novi_str1 += sequences[maxi]->data[ind1++];
            novi_str2 += sequences[maxj]->data[ind2++];
        }
        else if (pomstr[i] == 'U')
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