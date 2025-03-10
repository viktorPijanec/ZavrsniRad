#include <iostream>
#include "bioparser/fasta_parser.hpp"

using namespace std;

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

int Align(string seq1, string seq2)
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

    // cout << endl;
    // for (int i = 0; i < 10; ++i)
    // {
    //     for (int j = 0; j < 10; ++j)
    //     {
    //         cout << matrica[i][j].score << " ";
    //     }
    //     cout << '\n';
    // }

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
        cout << sequences[i]->name << endl;
        cout << sequences[i]->data << endl
             << endl;
    }

    for (int i = 0; i < sequences.size(); ++i)
    {
        for (int j = i + 1; j < sequences.size(); ++j)
        {
            cout << "Align(" << i + 1 << "," << j + 1 << ") : " << Align(sequences[i]->data, sequences[j]->data) << endl;
        }
    }
    return 0;
}