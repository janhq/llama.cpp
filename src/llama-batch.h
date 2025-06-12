#pragma once

#include "llama.h"

#include "llama-cparams.h"

#include <array>
#include <vector>
#include <set>
#include <bitset>
#include <unordered_map>

// keep this struct lightweight
// it points to data in `llama_batch_allocr`
struct llama_ubatch {
    bool equal_seqs;
    // TODO: whole_seqs for embeddings?

    uint32_t n_tokens;     // total tokens (n_seq_tokens * n_seqs)
    uint32_t n_seq_tokens; // tokens per sequence set
    uint32_t n_seqs;       // sequence sets in the ubatch
    uint32_t n_seqs_unq;   // unique sequence ids in the ubatch

    // seq_id_unq: unique sequence ids in the ubatch
    // seq_idx:    indices of the unique sequence ids in the ubatch in [0, n_seqs_unq)
    //             used for extracting sequence pooled embeddings

    //                          // size               | idx | val
    llama_token  *  token;      // [n_tokens]         | i   | id, token
    float        *  embd;       // [n_embd, n_tokens] | i   | embd
    llama_pos    *  pos;        // [n_tokens]         | i   | pos
    int32_t      *  n_seq_id;   // [n_tokens]         | i   | -
    llama_seq_id ** seq_id;     // [n_tokens]         | s   | s0, s1, seq_id
    llama_seq_id *  seq_id_unq; // [n_seqs_unq]       | s   | seq_id
    int32_t      *  seq_idx;    // [LLAMA_MAX_SEQ]    | -   | seq_idx
    int8_t       *  output;     // [n_tokens]         | i   | -
};

struct llama_sbatch_seq {
    int32_t n_seq_id;

    llama_seq_id * seq_id;

    size_t offset;
    size_t length;
};

// sequence-length-aware batch splitting
struct llama_sbatch {
    // tokens left in this batch
    size_t n_tokens;

    // only for debugging purposes
    const llama_vocab * vocab;

    // sorted indices into the batch
    std::vector<int64_t> ids;
    // batch indices of the output
    std::vector<int64_t> out_ids;
    std::vector<llama_sbatch_seq> seq;

    std::array<llama_seq_id, 1> seq_id_0 = { 0 }; // default sequence id

    // buffers for the ubatches
    // TODO: very hacky, this needs a complete rework
    struct ubatch_data {
        std::vector<llama_token>    token;
        std::vector<float>          embd;
        std::vector<llama_pos>      pos;
        std::vector<int32_t>        n_seq_id;
        std::vector<llama_seq_id *> seq_id;
        std::vector<int8_t>         output;
    };

    std::vector<ubatch_data> udatas;

    llama_ubatch reserve_ubatch(size_t n_ubatch, bool has_embd = false);

    void add_seq_to_ubatch(llama_ubatch & ubatch, llama_sbatch_seq & seq, size_t length);

    // simple split, unknown number of sequences of unequal lengths
    llama_ubatch split_simple(size_t n_ubatch);

    // make batches of equal-length sequences
    llama_ubatch split_equal(size_t n_ubatch);

    // sequence-wise split
    llama_ubatch split_seq(size_t n_ubatch);

    llama_sbatch() = default;
    llama_sbatch(const llama_batch & batch, size_t n_embd, bool simple_split = false);
};

    // used[i] indicates if token i has already been used in a previous ubatch
    std::vector<bool> used;

    std::array<llama_seq_id, 1> seq_id_0 = { 0 }; // default sequence id
    std::vector<llama_pos>      pos;
    std::vector<int32_t>        n_seq_id;
    std::vector<llama_seq_id *> seq_id;
    std::vector<int8_t>         output;

    int debug;
};
