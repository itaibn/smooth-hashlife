#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// I'm planning on using bignums for coordinates to allow sizes >2^64. Note that
// the current codebase does not support anything >2^31.
#include <gmp.h>

// GMP doesn't have a maximum or minimum function for mpz_t, and they'll be
// needed later
void my_mpz_max(mpz_t rop, const mpz_t op0, const mpz_t op1) {
    if (mpz_cmp(op0, op1) > 0) {
        mpz_set(rop, op0);
    } else {
        mpz_set(rop, op1);
    }
}

void my_mpz_min(mpz_t rop, const mpz_t op0, const mpz_t op1) {
    if (mpz_cmp(op0, op1) > 0) {
        mpz_set(rop, op1);
    } else {
        mpz_set(rop, op0);
    }
}

//// DEFINITIONS OF TYPES

struct block_struct;
typedef struct block_struct block;

// Leaf node for 4x4^H^H^H2x2 block (let's worry about efficiency later).
// Don't change this macro; this value is hardcoded in other places.
#define LEAFSIZE 2
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
} subblock;

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
    mp_bitcnt_t depth;
//    unsigned long refcount;
    block *res;
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

int
init_hash_cache() {
    int i, j, x, y;
    int p, tmp;
    unsigned long yhash, xyhash, hash;

    mpz_init_set_ui(hashprime_mpz, hashprime);

    yhash = 1;
    for (y=0; y<LEAFSIZE; y++) {
        xyhash = yhash;
        for (x=0; x<LEAFSIZE; x++) {
            point_hash_value[y*LEAFSIZE + x] = xyhash;
            xyhash = xyhash * xmul % hashprime;
        }
        yhash = yhash * ymul % hashprime;
    }

    for (p=0; p < 1 << (LEAFSIZE*LEAFSIZE); p++) {
        hash = 0;
        for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
            if (p & (1 << i)) {
                hash = (hash + 2*point_hash_value[i]) % hashprime;
            } else {
                hash = (hash + 1*point_hash_value[i]) % hashprime;
            }
        }
        leaf_hash_cache[p] = hash;
    }

    int size;
    unsigned long *table, point_value;
    for (x=1; x<=LEAFSIZE; x++) {
    for (y=1; y<=LEAFSIZE; y++) {
        size = 1 << (x*y);
        table = (unsigned long *) malloc(size * sizeof(unsigned long));
        for (p=0; p < size; p++) {
            hash = 0;
            for (i=0; i<x; i++) {
            for (j=0; j<y; j++) {
                point_value = point_hash_value[i + LEAFSIZE*j];
                if (p & (1 << (i+x*j))) {
                    hash = (hash + 2*point_value) % hashprime;
                } else {
                    hash = (hash + 1*point_value) % hashprime;
                }
                tmp = tmp >> 1;
            }}
            assert(hash < hashprime);
            table[p] = hash;
        }
        rectangle_hash_cache_size[(x-1) + LEAFSIZE*(y-1)] = size;
        rectangle_hash_cache[(x-1) + LEAFSIZE*(y-1)] = table;
    }
    }

    hash = point_hash_value[LEAFSIZE-1] * xmul % hashprime;
    for (i=0; i < 256; i++) {
        xmul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }

    hash = point_hash_value[LEAFSIZE*(LEAFSIZE-1)] * ymul % hashprime;
    for (i=0; i<256; i++) {
        ymul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }
}

// Note: d is the depth of the subnodes, not the combined node.
unsigned long
hash_node(unsigned long hnw, unsigned long hne, unsigned long hsw, unsigned long
        hse, mp_bitcnt_t d) {
    uint64_t xmul_d, ymul_d, xymul_d;
    xmul_d = xmul_cache[d];
    ymul_d = ymul_cache[d];
    xymul_d = xmul_d * ymul_d % hashprime;
    // Since hashprime is less than 2^31, all the terms in the sum are less than
    // 2^62, guaranteeing that there is no overflow.
    return (hnw + xmul_d*hne + ymul_d*hsw + xymul_d*hse) % hashprime;
}

// Here comes the most complicated function among the hashing utilities.
// Compute the hash of a rectangular subblock with northwest corner (x0, y0) and
// southeast corner (x1, y1). The rectangle is truncated if either (x0, y0) or
// (x1, y1) extend past the ends of base. If the parameter 'adjust' is set
// nonzero, the hash is calculated as if (x0, y0) is the origin.
unsigned long
hash_rectangle(block *base, const mpz_t ix0, const mpz_t ix1, const mpz_t iy0,
        const mpz_t iy1, int adjust) {
    //assert (mpz_cmp(ix0, ix1) =< 0 && mpz_cmp(iy0, y1) =< 0)

    unsigned long hash;
    mpz_t blocksize, zero, x0, x1, y0, y1;

    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(blocksize, LEAFSIZE);
    mpz_mul_2exp(blocksize, blocksize, base->depth);
    mpz_inits(x0, x1, y0, y1, NULL);
    my_mpz_max(x0, ix0, zero);
    my_mpz_max(y0, iy0, zero);
    my_mpz_min(x1, ix1, blocksize);
    my_mpz_min(y1, iy1, blocksize);
    //assert (mpz_sgn(x0) >= 0 && mpz_cmp(x1, blocksize) =< 0
    //    && mpz_sgn(y0) >= 0 && mpz_cmp(y1, blocksize) =< 0);
    
    if (mpz_cmp(x0, x1) >= 0 || mpz_cmp(y0, y1) >= 0) {
        hash = 0;
        goto end;
    }
    if (mpz_sgn(x0) == 0 && mpz_cmp(x1, blocksize) == 0 && mpz_sgn(y0) == 0 &&
        mpz_cmp(y1, blocksize) == 0) {
        hash = base->hash;
        goto end;
    }

    if (base->tag == LEAF_B) {
        unsigned long x0l, x1l, y0l, y1l;
        size_t size;
        unsigned long *table;
        leaf xmask, mask, pos, row, rect;
        int i;
        x0l = mpz_get_ui(x0);
        x1l = mpz_get_ui(x1);
        y0l = mpz_get_ui(y0);
        y1l = mpz_get_ui(y1);
        i = (x1l-x0l-1) + LEAFSIZE*(y1l-y0l-1);
        assert (0 <= i && i < LEAFSIZE*LEAFSIZE);
        table = rectangle_hash_cache[i];
        size = rectangle_hash_cache_size[i];
        pos = base->content.b_l;
        rect = 0;
        xmask = ((1 << (x1l - x0l)) - 1) << x0l;
        for (i = y0l; i < y1l; i++) {
            mask = xmask << (i*LEAFSIZE);
            row = (pos & mask) >> (i*LEAFSIZE+x0l);
            rect |= row << ((i-y0l)*(x1l-x0l));
        }
        assert(rect < 1<<((x1l-x0l)*(y1l-y0l)));
        assert(rect < size);
        assert(size == 1<<((x1l-x0l)*(y1l-y0l)));
        hash = point_hash_value[x0l + LEAFSIZE*y0l] * table[rect] % hashprime;
        printf(" h %lu\n", hash);
        goto end;
    } else if (base->tag != NODE_B) {
        fprintf(stderr, "CONTAIN_B not supported as input to hash_rectangle()");
    }

    node n = base->content.b_n;
    unsigned long hnw, hne, hsw, hse; 
    mpz_t halfblock;

    mpz_t shiftx0, shiftx1, shifty0, shifty1;
    mpz_inits(halfblock, shiftx0, shiftx1, shifty0, shifty1, NULL);
    mpz_tdiv_q_2exp(halfblock, blocksize, 1);
    mpz_sub(shiftx0, x0, halfblock);
    mpz_sub(shiftx1, x1, halfblock);
    mpz_sub(shifty0, y0, halfblock);
    mpz_sub(shifty1, y1, halfblock);
    hnw = hash_rectangle(NW(n), x0, x1, y0, y1, 0);
    hne = hash_rectangle(NE(n), shiftx0, shiftx1, y0, y1, 0);
    hsw = hash_rectangle(SW(n), x0, x1, shifty0, shifty1, 0);
    hse = hash_rectangle(SE(n), shiftx0, shiftx1, shifty0, shifty1, 0);
    
    hash = hash_node(hnw, hne, hsw, hse, base->depth-1);
    mpz_clears(halfblock, shiftx0, shiftx1, shifty0, shifty1, NULL);

    printf("[%lu %lu %lu %lu -> %lu]\n", hnw, hne, hsw, hse, hash);
end:
    if (adjust) {
        mpz_t x_adj, y_adj;
        uint64_t xy_adj;
        mpz_init_set_ui(x_adj, xmul);
        mpz_neg(x0, x0);
        mpz_powm(x_adj, x_adj, x0, hashprime_mpz);
        mpz_init_set_ui(y_adj, ymul);
        mpz_neg(y0, y0);
        mpz_powm(y_adj, y_adj, y0, hashprime_mpz);
        xy_adj = mpz_get_ui(x_adj) * mpz_get_ui(y_adj) % hashprime;
        hash = (unsigned long) ((xy_adj * (uint64_t) hash) % hashprime);
        mpz_clears(x_adj, y_adj, NULL);
    }
    mpz_clears(blocksize, zero, x0, x1, y0, y1, NULL);
    return hash;
}

//// BASIC BLOCK CREATION FUNCTIONS

// Allocate a block with a given hash from the hash table
block * new_block(unsigned long hash);

block *
mkblock_leaf(leaf l) {
    unsigned long hash = leaf_hash_cache[l];
    block *b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = LEAF_B;
    b->content.b_l = l;
    b->depth = 0;
    return b;
}

// Combine four blocks into a node block
block *
mkblock_node(block *nw, block *ne, block *sw, block *se) {
    unsigned long hash;
    block *b;
    mp_bitcnt_t d;

    if (nw == NULL || ne == NULL || sw == NULL || se == NULL) {
        return NULL;
    }
    d = nw->depth;
    if (ne->depth != d || sw->depth != d || se->depth != d) {
        return NULL;
    }
    if (d >= 256) {
        fprintf(stderr, "This implementation currently does not supported sizes"
            "larger than 2^257\n");
        return NULL;
    }

    node n = {{{nw, ne}, {sw, se}}};
    hash = hash_node(nw->hash, ne->hash, sw->hash, se->hash, d);
    b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = NODE_B;
    b->content.b_n = n;
    b->depth = d+1;
    return b;
}

block *
mkblock_contain(block *superblock, mpz_t x, mpz_t y) {
    mp_bitcnt_t d;
    unsigned long hash;
    mpz_t size, x_east, y_south;
    block *b;

    d = superblock->depth-1;
    assert(d >= 0);
    mpz_init_set_ui(size, LEAFSIZE);
    mpz_mul_2exp(size, size, d);
    mpz_init(x_east);
    mpz_add(x_east, x, size);
    mpz_init(y_south);
    mpz_add(y_south, y, size);

    assert(mpz_sgn(x) >= 0 && mpz_cmp(x_east, size) <= 0 &&
        mpz_sgn(y) >= 0 && mpz_cmp(y_south, size) <= 0
        && superblock);

    hash = hash_rectangle(superblock, x, x_east, y, y_south, 1);

    b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = CONTAIN_B;
    b->depth = d;
    b->content.b_c.superblock = superblock;
    mpz_init_set(b->content.b_c.x, x);
    mpz_init_set(b->content.b_c.y, y);
    // ... b->content.b_c.t ...
    mpz_clears(size, x_east, y_south, NULL);
    return b;
}

// Given a 2^nx2^n block b, return the 2^(n-1)x2^(n-1) subblock with northwest
// corner at (i*2^(n-2), j*2^(n-2)). For example:
//
//      N
//  +-+---+-+
//  | |   | |
//  | | r | |
//  | |   | |
// W| +---+ |E    r = block_index(A, 1, 0)
//  |   A   |
//  |       |
//  |       |
//  +-------+
//      S
block *
block_index(block *b, int i, int j) {
    if (b == NULL) {return NULL;}
    if (b->tag == LEAF_B) {return NULL;}

    if (b->tag == CONTAIN_B) {
        mpz_t x, y, halfblock, shift;
        mpz_inits(x, y, halfblock, shift, NULL);
        // TODO
        mpz_clears(x, y, halfblock, shift, NULL);
    } else if (b->tag != NODE_B) {
        fprintf(stderr, "Invalid tag %d for block at %p (hash %lu)\n", b->tag,
            b, b->hash);
    }

    node b_no = b->content.b_n;
    
    if (((i | j) & 1) == 0) {
        return CORNER(b_no, i, j);
    } else {
        node n, tmpn;
        block *tmpb;
        int p, i0, i1, q, j0, j1;
        for (p=0; p<2; p++) {
        for (q=0; q<2; q++) {
            i0 = ((p + i) & 2) >> 1;
            i1 = (p + i) & 1;
            j0 = ((q + j) & 2) >> 1;
            j1 = (q + j) & 1;
            tmpb = CORNER(b_no, i0, j0);
            tmpn = tmpb->content.b_n;
            CORNER(n, p, q) = CORNER(tmpn, i1, j1);
        }}
        return mkblock_node(NW(n), NE(n), SW(n), SE(n));
    }
}

//// CA COMPUTATION PROPER

// Lookup table for 1-step of 4x4 blocks. Idea stolen from Tomas Rokicki's
// hlife.c.
leaf result[65536];

void
init_result() {
    int result33[65536];
    int pos, i, bitcnt, tmp;
    for (i=0; i<512; i++) {
        pos = (i&0b111) | ((i&0b111000) << 1) | ((i&0b111000000) << 2);
        tmp = i & ~0x20;
        while (tmp > 0) {
            bitcnt += tmp & 1;
            tmp = tmp >> 1;
        }
        bitcnt = __builtin_popcount(pos & ~0x20);
        int res;
        if (i & 0x10) {
            res = bitcnt == 2 || bitcnt == 3;
        } else {
            res = bitcnt == 3;
        }
        if (res) {
            result33[pos] = 1;
        } else {
            result33[pos] = 0;
        }
    }
    for (i=0; i<65536; i++) {
        result[i] = result33[i&0x0777] |
            (result33[(i&0x0eee) >> 1] << 1) |
            (result33[(i&0x7770) >> 4] << 2) |
            (result33[(i&0xeee0) >> 5] << 3);
    }
}

// Sorry, the code here is a bit messy and repetitive.
block *
evolve(block *x) {
    block *r;
    if (x == NULL) {
        fprintf(stderr, "NULL input to evolve()\n");
    }
    if (x->res) {
        return x->res;
    }
    if (x->tag != NODE_B) {
        fprintf(stderr, "Only nodes have evolve() implemented\n");
        return NULL;
    }
    if (x->depth == 1) {
        // Subblocks are all leafs
        node n = x->content.b_n;
        assert(NW(n)->tag == LEAF_B);
        int unpack_x = NW(n)->content.b_l & 3 |
                       (NE(n)->content.b_l & 3) << 2 |
                       (NW(n)->content.b_l & 12) << 2 |
                       (NE(n)->content.b_l & 12) << 4 |
                       (SW(n)->content.b_l & 3) << 8 |
                       (SE(n)->content.b_l & 3) << 10 |
                       (SW(n)->content.b_l & 12) << 10 |
                       (SE(n)->content.b_l & 12) << 12;
        r = mkblock_leaf(result[unpack_x]);
    } else {
        // Half-sized subblocks of x on the north, south, east, west, and
        // center:
        block *n, *s, *w, *e, *c;
        node no = x->content.b_n;
        assert(NW(no)->tag == NODE_B);
        node tmp; // Not an actual node; just a convient way to store four
                  // blocks associated with corners.
        // Recall mkblock_node(nw, ne, sw, se)
        // This part is tedious and error-prone
        n = mkblock_node(NE(NW(no)->content.b_n),
                         NW(NE(no)->content.b_n),
                         SE(NW(no)->content.b_n),
                         SW(NE(no)->content.b_n));
        s = mkblock_node(NE(SW(no)->content.b_n),
                         NW(SE(no)->content.b_n),
                         SE(SW(no)->content.b_n),
                         SW(SE(no)->content.b_n));
        w = mkblock_node(SW(NW(no)->content.b_n),
                         SE(NW(no)->content.b_n),
                         NW(SW(no)->content.b_n),
                         NE(SW(no)->content.b_n));
        e = mkblock_node(SW(NE(no)->content.b_n),
                         SE(NE(no)->content.b_n),
                         NW(SE(no)->content.b_n),
                         NE(SE(no)->content.b_n));
        c = mkblock_node(SE(NW(no)->content.b_n),
                         SW(NE(no)->content.b_n),
                         NE(SW(no)->content.b_n),
                         NW(SE(no)->content.b_n));
        n = evolve(n);
        s = evolve(s);
        w = evolve(w);
        e = evolve(e);
        c = evolve(c);
        NW(tmp) = evolve(mkblock_node(evolve(NW(no)), n, w, c));
        NE(tmp) = evolve(mkblock_node(n, evolve(NE(no)), c, e));
        SW(tmp) = evolve(mkblock_node(w, c, evolve(SW(no)), s));
        SE(tmp) = evolve(mkblock_node(c, e, s, evolve(SE(no))));
        r = mkblock_node(NW(tmp), NE(tmp), SW(tmp), SE(tmp));
    }
    return (x->res = r);
}


//// READING AND WRITING BLOCKS

block *write_bit(block *b, unsigned long y, unsigned long x, char bit);

block *
read_life_105(FILE *f) {
    int c;
    char bit;

    while ((c = fgetc(f)) == '#') {
        if (c == EOF) return NULL;
        while ((c = fgetc(f)) != '\n') {
            if (c == EOF) return NULL;
        }
    }
    
    unsigned long y, x;
    int res;
    block *b, *empty, *tmp;
    b = empty = mkblock_leaf(0);
    x = y = 0;
    while(1) {
        switch (c) {
            case '*':
                bit = 1;
                goto star_or_dot;
            case '.':
                bit = 0;
            star_or_dot:
                while (!(tmp = write_bit(b, y, x, bit))) {
                    b = mkblock_node(b, empty, empty, empty);
                    empty = mkblock_node(empty, empty, empty, empty);
                }
                b = tmp;
                x++;
                break;
            case '\r':
                if ((c = fgetc(f)) != '\n') {
                    return NULL;
                }
            case '\n':
                x = 0;
                y++;
                break;
            case EOF:
                return b;
                break;
            default:
                return NULL;
                break;
        }
        c = fgetc(f);
    }
}

block *
read_mc(FILE *f) {
    return NULL; // Screw this. I should probably only start working on reading
                 // .mc files after I have the hash-table code to support it.

    char c;

    // Read first line: "[M2] ..."
    while ((c = getc(f)) != '\n') {
        if (c == EOF) {
            return NULL;
        }
    } 

    block *blank_eight, **blocktable;
    blank_eight = mkblock_leaf(0);
    blank_eight = mkblock_node(blank_eight, blank_eight, blank_eight,
        blank_eight);
    blank_eight = mkblock_node(blank_eight, blank_eight, blank_eight,
        blank_eight);
    if (LEAFSIZE != 2)
        {fprintf(stderr, "read_mc() needs LEAFSIZE=2\n"); return NULL;}
    blocktable /*= ...*/;

    while(1) {
        c = getc(f);
        if (c == '$' || c == '*' || c == '.') {
            int x=0, y=0;
            while (c != '\n') {
                switch(c) {
                    //case '
                    }}}}
}

// Used in read_life_105. Returns NULL when (x,y) is out of range.
block *
write_bit(block *b, unsigned long y, unsigned long x, char bit) {
    unsigned long size = 2;
    /*
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    */
    size = 2 << b->depth;

    if (x >= size) {
        //return -1;
        return NULL;
    }
    if (y >= size) {
        //return -2;
        return NULL;
    }

    if (b->tag == LEAF_B) {
        uint16_t mask = 1 << (2*y+x);
        leaf new_leaf = (b->content.b_l & ~mask) | (bit ? mask : 0);
        return mkblock_leaf(new_leaf);
    } else if (b->tag == NODE_B) {
        node n = b->content.b_n;
        if (y < size/2) {
            if (x < size/2) {
                NW(n) = write_bit(NW(n), y, x, bit);
            } else {
                NE(n) = write_bit(NE(n), y, x-size/2, bit);
            }
        } else {
            if (x < size/2) {
                SW(n) = write_bit(SW(n), y-size/2, x, bit);
            } else {
                SE(n) = write_bit(SE(n), y-size/2, x-size/2, bit);
            }
        }
        if (NW(n) == NULL || NE(n) == NULL || SW(n) == NULL || SE(n) == NULL) {
            return NULL;
        } else {
            return mkblock_node(NW(n), NE(n), SW(n), SE(n));
        }
    } else {
        fprintf(stderr, "No code for writing in a block neither leaf nor node");
        return NULL;
    }
}

int
print_line(block *b, long y, FILE *f) {
    unsigned long size = 2;
    /*
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    */
    size = 2 << b->depth;

    if (b->tag == LEAF_B) {
        int out = 0xf & b->content.b_l >> (2*y);
        int j;
        for (j = 0; j<2; j++) {
            fputc((out&1) ? '*' : '.', f);
            out >>= 1;
        }
    } else if (b->tag == NODE_B) {
        if (y < size/2) {
            print_line(NW(b->content.b_n), y, f);
            print_line(NE(b->content.b_n), y, f);
        } else {
            print_line(SW(b->content.b_n), y - size/2, f);
            print_line(SE(b->content.b_n), y - size/2, f);
        }
    }
}

int
display(block *b, FILE *f) {
    unsigned long size = 2;
    /*
    block *tmp = b;
    while (tmp->tag != LEAF_B) {
        assert(tmp->tag == NODE_B);
        size <<= 1;
        tmp = tmp->content.b_n.nw;
    }
    */
    size = 2 << b->depth;
    if (size == 0) {
        fprintf(stderr, "Overflow error: code doesn't support blocks of this"
            "size\n");
        return 1;
    }

    long i;
    for (i = 0; i < size; i++) {
        print_line(b, i, f);
        fputc('\n', f);
    }
}

//// HASH TABLE-RELATED FUNCTIONS
// Note: Currently this sucks, there is no hash table resizing nor garbage
// collection.

unsigned int ht_size = 1000000;
block **hashtable;

void
init_hashtable() {
    hashtable = (block **) malloc(ht_size * sizeof(block *));
    memset(hashtable, 0, ht_size * sizeof(block *));
}

block *
new_block(unsigned long hash) {
    unsigned long ii = 0;
    int i;
    block *b;
    do {
        ii = ii + hash;
        i = ii % ht_size;
        b = hashtable[i];
        if (b == NULL) {
            b = (block *) malloc(sizeof(block));
            b->hash = hash;
            b->res = NULL;
            b->tag = EMPTY;
//            b->refcount = 0;
            // Insert command here of the form:
            //  num_entries++;
        }
    } while (b->hash != hash);
//    b->refcount++;
    return b;
}

/*
// Decrease the refrence count of a block.
int
unref(block *b) {
    if (--(b->refcount)) {
        return 1;
    } else {
        
*/

//// MAIN

main(int argc, char **argv) {
    init_hashtable();
    init_hash_cache();
    init_result();

    if (argc != 6) {
        fprintf(stderr, "shlife pattern x0 x1 y0 y1\n");
        exit(1);
    }

    FILE *f;
    block *b;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(1);
    }
    b = read_life_105(f);
    if (b == NULL) {
        fprintf(stderr, "Badly formatted input\n");
        exit(1);
    }
    fclose(f);

    mpz_t x[4];
    int i, err=0;
    for (i=0; i<4; i++) {
        err |= mpz_init_set_str(x[i], argv[i+2], 10);
    }
    if (err) {
        fprintf(stderr, "Arguments x0, x1, y0, y1 must be numbers\n");
        exit(1);
    }

    //display(evolve(b), stdout);
    display(b, stdout);
    //if (b->tag == NODE_B) display(b->content.b_n.nw, stdout);
    printf("%lu\n", b->hash);
    printf("%lu\n", hash_rectangle(b, x[0], x[1], x[2], x[3], 1));
    exit(0);
}
