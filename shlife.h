#ifndef H_SHLIFE
#define H_SHLIFE

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// I'm planning on using bignums for coordinates to allow sizes >2^64. Note that
// the current codebase does not support anything >2^31.
#include <gmp.h>
void my_mpz_min(mpz_t rop, const mpz_t op0, const mpz_t op1);
void my_mpz_max(mpz_t rop, const mpz_t op0, const mpz_t op1);

#ifdef DEBUG
#define TRACE gmp_printf
#define DEBUG 1
#else
#define TRACE(a, ...)
#define DEBUG 0
#endif // DEBUG

//// DEFINITIONS OF TYPES

struct block_struct;
typedef struct block_struct block;

block **hashtable;

// Integer type for the depth of the block-tree.
typedef mp_bitcnt_t depth_t;

// Leaf node for 4x4^H^H^H2x2 block (let's worry about efficiency later).
// Don't change this macro; this value is hardcoded in other places.
#define LEAFSIZE 2
#define LGLEAFSIZE 1
/*
Data format is little-endian going west-to-east, north-to-south, like so:
   N

  0 1
W     E
  2 3

   S
*/
typedef uint8_t leaf;

// A node block recursively containing 4 smaller blocks.
//typedef block *(node[2][2]);
typedef struct {
    block *corner[2][2];
} node;

#define NW(n) ((n).corner[0][0])
#define NE(n) ((n).corner[0][1])
#define SW(n) ((n).corner[1][0])
#define SE(n) ((n).corner[1][1])
#define CORNER(n, i, j) ((n).corner[i][j])

typedef struct {
    mpz_t x, y/*, t*/;
    block *superblock;
    block *as_node;
} subblock;

struct inner_pattern;

// General block. May be a leaf or node.
struct block_struct {
    enum b_tag {
        EMPTY = 0, // For uninitialized block entries in the hash table.
        LEAF_B,
        NODE_B,
        CONTAIN_B
    } tag;
    union {
        leaf b_l;
        node b_n;
        subblock b_c;
    } content;
    unsigned long hash;
    depth_t depth;
//    unsigned long refcount;
    block *res;
    
    int nfocus;
    struct inner_pattern *foci;

    int index;
};

#define LGLENGTH(b) ((b)->depth+LGLEAFSIZE)

struct inner_pattern {
    mpz_t x, y /*, t*/;
    block *pattern;
    int depth_diff;
};

//// TOOLS FOR CALCULATING THE HASH FUNCTION
// The intended hash function is as follows: Given a block b of size 2^nx2^n, it
// sums modulo hashprime for all 0 =< i,j < 2^n the value xmul^i*ymul^j if the
// coordinate (i,j) is dead in b and 2*xmul^i*ymul^j if (i, j) is alive in b.
// This way simply multiplying the hash by xmul^n*ymul^m "translates" the block
// by (n, m), which allows calculating the hashes of subblocks more easily.

// Currently very little thought was put in choosing hashprime, xmul, & ymul.
// hashprime was chosen as a large prime that fits in 32 bits, and xmul & ymul
// were chosen so as not to have any obvious linear relationships which might
// cause collisions. It would be good to have that both xmul and ymul are
// primitive roots but this hasn't been checked.

// Note: hashprime must be less than 2^31 for code to work. See comment in
// hash_node().

const unsigned long hashprime = 1000000007;
mpz_t hashprime_mpz;
const unsigned long xmul = 2331, ymul = 121212121;

// {x,y}mul_cache[i] caches the value {xmul,ymul}^(LEAFSIZE*2^i), and is used
// for computing the hash of a depth i node. Currently no support for depths
// >=256.
unsigned long xmul_cache[256];
unsigned long ymul_cache[256];

// A cache of the hashes of all leaf nodes.
unsigned long leaf_hash_cache[1 << (LEAFSIZE*LEAFSIZE)];
// A set of caches for hashes of rectangular subblocks of leaf nodes of all
// sizes.
size_t rectangle_hash_cache_size[LEAFSIZE * LEAFSIZE];
unsigned long *rectangle_hash_cache[LEAFSIZE * LEAFSIZE];
// Values xmul^i*ymul^j for 0 =< i,j < LEAFSIZE
unsigned long point_hash_value[LEAFSIZE * LEAFSIZE];

// Initialize all the caches defined above
int init_hash_cache();

// Given the hashes of all the corners of a node, and the depth of these
// corners, calculate the hash of the whole node.
unsigned long hash_node(unsigned long hnw, unsigned long hne, unsigned long hsw,
    unsigned long hse, depth_t d);

// Here comes the most complicated function among the hashing utilities.
// Compute the hash of a rectangular subblock with northwest corner (x0, y0) and
// southeast corner (x1, y1). The rectangle is truncated if either (x0, y0) or
// (x1, y1) extend past the ends of base. If the parameter 'adjust' is set
// nonzero, the hash is calculated as if (x0, y0) is the origin.
unsigned long hash_rectangle(block *base, const mpz_t ix0, const mpz_t ix1,
    const mpz_t iy0, const mpz_t iy1, int adjust);

//// HASH TABLE-RELATED FUNCTIONS
// Note: Currently this sucks, there is no hash table resizing nor garbage
// collection.

unsigned int ht_size = 1000000;
block **hashtable;

int init_hashtable();

block *new_block(unsigned long hash);

/*
// Decrease the refrence count of a block.
int unref(block *b);
*/

//// BASIC BLOCK CREATION FUNCTIONS

// Allocate a block with a given hash from the hash table
block *new_block(unsigned long hash);

// Make a block out of a leaf
block *mkblock_leaf(leaf l);

// Combine four blocks into a node block
block *mkblock_node(block *nw, block *ne, block *sw, block *se);

// TODO: Seperate this into two different functions, one for diff=1, one for
// diff>1 (possibly on diff=2 necessary). Thus far this procedure is a mixup
// between block creation and more complicated block processing.
block *mkblock_contain(block *superblock, mpz_t x, mpz_t y, depth_t diff, int
    rec);

// Given a 2^nx2^n block b, return the 2^(n-1)x2^(n-1) subblock with northwest
// corner at (jx*2^(n-2), iy*2^(n-2)). For example:
//
//      N
//  +-+---+-+
//  | |   | |
//  | | r | |
//  | |   | |
// W| +---+ |E    r = block_index(A, 0, 1)
//  |   A   |
//  |       |
//  |       |
//  +-------+
//      S
block *
block_index(block *b, int iy, int jx);

// Make a block with every cell dead of sidelength 2^lglength
block *blank_block(depth_t lglength);

// Write over coordinate (x,y) of b to have liveness bit, and return the
// resulting block. (This does not actually modify the block b.)
block *write_bit(block *b, unsigned long y, unsigned long x, char bit);

//// CA COMPUTATION PROPER

int add_foci();

// Lookup table for 1-step of 4x4 blocks. Idea stolen from Tomas Rokicki's
// hlife.c.
leaf result[65536];

int init_result();

block *evolve(block *x);

block *half_evolve(block *b, int i, int j);

//// READING AND WRITING BLOCKS

block *read_life_105(FILE *f);

int display_raw(block *b, FILE *f);

block *read_mc(FILE *f);

// Useful for debugging.
int display_raw(block *b, FILE *f);

// Used in read_life_105. Returns NULL when (x,y) is out of range.
block *write_bit(block *b, unsigned long y, unsigned long x, char bit);

// Display the line from (x0, y) to (x1, y) in Life 1.05 format onto file f.
// Ends with a newline if newline!=0.
int display_line(block *b, mpz_t y, mpz_t x0, mpz_t x1, int newline, FILE *f);

// Displays a block in Life 1.05 format.
int display(block *b, FILE *f);

#endif //H_SHLIFE
