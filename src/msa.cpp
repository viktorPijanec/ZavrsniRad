#include <iostream>
#include <unordered_map>
#include "bioparser/fasta_parser.hpp"

using namespace std;

namespace std
{
    template <>
    struct hash<std::pair<int, int>>
    {
        std::size_t operator()(const std::pair<int, int> &p) const noexcept
        {
            return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
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

    unordered_map<pair<int, int>, int> alignments;

    for (const auto &[key, value] : alignments)
    {
        cout << key.first << "," << key.second << " : " << value << endl;
    }
    int maxi = 0;
    int maxj = 0;
    int maxalign = -100000;
    for (int i = 0; i < sequences.size(); ++i)
    {
        for (int j = i + 1; j < sequences.size(); ++j)
        {
            pair<int, int> par = make_pair(i, j);
            alignments.insert({par, Align(sequences[i]->data, sequences[j]->data, nullptr)});
            if (alignments[par] > maxalign)
            {
                maxi = i;
                maxj = j;
                maxalign = alignments[par];
            }
            cout << "Align(" << i + 1 << "," << j + 1 << ") : " << alignments[par] << endl;
        }
    }
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