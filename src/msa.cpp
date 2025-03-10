#include <iostream>
#include "bioparser/fasta_parser.hpp"

using namespace std;

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

int Align(string seq1, string seq2) {
    
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

    return 0;
}