#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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
typedef struct {
    block *nw, *ne, *sw, *se;
} node;

typedef struct {
    mpz_t x, y/*, t*/;
    block *superblock;
} subblock;

// General block. May be a leaf or node.
// To be defined: subblocks of bigger blocks.
struct block_struct {
    enum b_tag {
        //EMPTY, // For uninitialized block entries in the hash table.
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
// mkblock_node.
const unsigned long hashprime = 1000000007;
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
unsigned long *rectangle_hash_cache[(LEAFSIZE-1)*(LEAFSIZE-1)];
// Values xmul^i*ymul^j for 0 =< i,j < LEAFSIZE
unsigned long point_hash_value[LEAFSIZE * LEAFSIZE];

int
init_hash_cache() {
    int i, x, y, hash;
    int p, tmp;
    unsigned long yhash, xyhash;
    yhash = 1;
    for (x=0; x<LEAFSIZE; x++) {
        xyhash = yhash;
        for (y=0; y<LEAFSIZE; y++) {
            printf("%d: %lu\n", x*LEAFSIZE + y, xyhash);
            point_hash_value[x*LEAFSIZE + y] = xyhash;
            xyhash = xyhash * xmul % hashprime;
        }
        yhash = yhash * ymul % hashprime;
    }
    for (p=0; p < 1 << (LEAFSIZE*LEAFSIZE); p++) {
        hash = 0;
        tmp = p;
        for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
            if (tmp & 1) {
                hash = (hash + 2*point_hash_value[i]) % hashprime;
            } else {
                hash = (hash + 1*point_hash_value[i]) % hashprime;
            }
            tmp = tmp >> 1;
        }
        leaf_hash_cache[p] = hash;
    }

    int size;
    unsigned long *table;
    for (x=1; x<LEAFSIZE; x++) {
    for (y=1; y<LEAFSIZE; y++) {
        size = 1 << (x*y);
        table = (unsigned long *) malloc(size * sizeof(unsigned long));
        for (p=0; p < 1 << (LEAFSIZE*LEAFSIZE); p++) {
            hash = 0;
            tmp = p;
            for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
                if (tmp & 1) {
                    hash = (hash + 2*point_hash_value[i]) % hashprime;
                } else {
                    hash = (hash + 1*point_hash_value[i]) % hashprime;
                }
                tmp = tmp >> 1;
            }
        }
        rectangle_hash_cache[(x-1) + (LEAFSIZE-1)*(y-1)] = table;
    }
    }

    hash = point_hash_value[LEAFSIZE-1] * xmul % hashprime;
    printf("x: %lu\n", hash);
    for (i=0; i < 256; i++) {
        xmul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }

    hash = point_hash_value[LEAFSIZE*(LEAFSIZE-1)] * ymul % hashprime;
    printf("y: %lu\n", hash);
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

/*
// Meaning of corner: 0->nw; 1->sw; 2->ne; 3->se;
unsigned long
corner_hash(block *base, mpz_t x, mpz_t y, int corner) {
    if (base->depth == 0) {
        // ???
    }

    if (base->tag == CONTAIN_B) {
        fprintf(stderr, "no support for contain type blocks in"
            "'corner_hash'\n");
        return hashprime;
    }

    node n = base->content.b_n
    block *ii, *io, *oi, *oo;

    switch (corner) {
        case 0:
            ii = n.nw; io = n.ne; oi = n.sw; oo = n.se;
            break;
        case 1:
            ii = n.sw; io = n.se; oi = n.nw; oo = n.ne;
            break;
        case 2:
            ii = n.ne; io = n.nw; oi = n.se; oo = n.sw;
            break;
        case 3:
            ii = n.se; io = n.sw; oi = n.ne; oo = n.nw;
            break;
        default:
            fprintf(stderr, "Invalid parameter corner to 'corner_hash': %d\n",
                corner);
            return hashprime;
    }

    mpz_t halfblock;
    int cmp0, cmp1;
    mpz_init_set_ui(halfblock, LEAFSIZE);
    mpz_mul_2exp(halfblock, halfblock, b->depth - 1);
    cmp0 = 
*/

// Compute the hash of a rectangular subblock with northwest corner (x0, y0) and
// southeast corner (x1, y1). The rectangle is truncated if 
// XXX UNFINISHED
unsigned long
hash_rectangle(block *base, const mpz_t ix0, const mpz_t ix1, const mpz_t iy0,
        const mpz_t iy1) {
    //assert (mpz_cmp(ix0, ix1) =< 0 && mpz_cmp(iy0, y1) =< 0)

    unsigned long hash;
    mpz_t blocksize, zero, x0, x1, y0, y1;

    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(blocksize, LEAFSIZE);
    mpz_mul_2exp(blocksize, blocksize, base->depth);
    my_mpz_max(x0, ix0, zero);
    my_mpz_max(y0, iy0, zero);
    my_mpz_min(x1, ix1, blocksize);
    my_mpz_min(y1, iy1, blocksize);
    //assert (mpz_sgn(x0) >= 0 && mpz_cmp(x1, blocksize) =< 0
    //    && mpz_sgn(y0) >= 0 && mpz_cmp(y1, blocksize) =< 0);
    
    if (mpz_cmp(x0, x1) <= 0 || mpz_cmp(y0, y1) <= 0) {
        return 0;
    }
    if (mpz_sgn(x0) == 0 && mpz_cmp(x1, blocksize) == 0 && mpz_sgn(y0) == 0 &&
        mpz_cmp(y1, blocksize) == 0) {
        return base->hash;
    }

    if (base->tag == LEAF_B) {
        unsigned long x0l, x1l, y0l, y1l;
        unsigned long *table;
        leaf xmask, pos;
        int i;
        x0l = mpz_get_ui(x0);
        x1l = mpz_get_ui(x1);
        y0l = mpz_get_ui(y0);
        y1l = mpz_get_ui(y1);
        //hash = point_hash_value[x0l + LEAFSIZE*y0l]
        //    * rectangle_hash_cache[(x1l-x0l-1) + (LEAFSIZE-1)*(y1l-y0l-1)]
        //    % hashprime;
        table = rectangle_hash_cache[(x1l-x0l-1) + (LEAFSIZE-1)*(y1l-y0l-1)];
        xmask = (1 << (x1l - x0l) - 1) << x0l;
        pos = 0;
        for (i = y0l; i < y1l; i++) {
            pos |= xmask << ((x1l-x0l)*i);
        }
        hash = table[hash];
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
    hnw = hash_rectangle(n.nw, x0, x1, y0, y1);
    hne = hash_rectangle(n.ne, shiftx0, shiftx1, y0, y1);
    hsw = hash_rectangle(n.sw, x0, x1, shifty0, shifty1);
    hse = hash_rectangle(n.se, shiftx0, shiftx1, shifty0, shifty1);
    
    hash = hash_node(hnw, hne, hsw, hse, base->depth-1);
    mpz_clears(halfblock, shiftx0, shiftx1, shifty0, shifty1, NULL);
end:
    mpz_clears(blocksize, zero, x0, x1, y0, y1, NULL);
    return hash;

    /*
    mpz_t tmpx0, tmpx1, tmpy0, tmpy1;
    mpz_init(halfblock);
    mpz_init(tmpx0); mpz_init(tmpx1); mpz_init(tmpy0); mpz_init(tmpy1);
    mpz_tdiv_q_2exp(halfblock, blocksize, 1);

    if (mpz_cmp(x0, halfblock) < 0 && mpz_cmp(y0, halfblock) < 0) {
        my_mpz_min(tmpx, x1, halfblock);
        my_mpz_min(tmpy, y1, halfblock);
        hnw = hash_rectangle(n.nw, x0, tmpx, y0, tmpy);
    } else {
        hnw = 0;
    }
    // Copy & Paste code *Sigh*
    if (mpz_cmp(x1, halfblock) > 0 && mpz_cmp(y0, halfblock) < 0) {
        my_mpz_max(tmpx, x0, halfblock);
        my_mpz_min(tmpy, y1, halfblock);
        hne = hash_rectangle(n.ne, tmpx, x1, y0, tmpy);
    } else {
        hnw = 0;
    }
    if (mpz_cmp(x0, halfblock) < 0 && mpz_cmp(y1, halfblock) > 0) {
        my_mpz_min(tmpx, x1, halfblock);
        my_mpz_max(tmpy, y0, halfblock);
        hsw = hash_rectangle(n.
    } else {
        hnw = 0;
    }
    if (mpz_cmp(x1, halfblock) > 0 && mpz_cmp(y1, halfblock) > 0) {
        my_mpz_max(tmpx, x0, halfblock);
        my_mpz_max(tmpy, y0, halfblock);
        hnw = hash_rectangle(n.nw, x0, tmpx, y0, tmpy);
    } else {
        hnw = 0;
    }
    */
}

//// BASIC BLOCK CREATION FUNCTIONS

// Allocate a block with a given hash from the hash table
block * new_block(unsigned long hash);

block *
mkblock_leaf(leaf l) {
    // Hash function will be changed later
    //unsigned long hash = (l*l*l+94455) % hashprime;
    unsigned long hash = leaf_hash_cache[l];
    block *b = new_block(hash);
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
    //d++;

    node n = {nw, ne, sw, se};
    // Calculate the hash function
    // Moved to hash_node
    /*
    uint64_t xmul_d, ymul_d, xymul_d, nwh, neh, swh, seh;
    nwh = nw->hash; neh = ne->hash; swh = sw->hash; seh = se->hash;
    xmul_d = xmul_cache[d];
    ymul_d = ymul_cache[d];
    xymul_d = xmul_d * ymul_d % hashprime;
    // Since hashprime is less than 2^31, all the terms in the sum are less than
    // 2^62, guaranteeing that there is no overflow.
    hash = (nwh + xmul_d*neh + ymul_d*swh + xymul_d*seh) % hashprime;
    */
    hash = hash_node(nw->hash, ne->hash, sw->hash, se->hash, d);
    b = new_block(hash);
    b->tag = NODE_B;
    b->content.b_n = n;
    b->depth = d+1;
    return b;
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
        assert(n.nw->tag == LEAF_B);
        int unpack_x = n.nw->content.b_l & 3 |
                       (n.ne->content.b_l & 3) << 2 |
                       (n.nw->content.b_l & 12) << 2 |
                       (n.ne->content.b_l & 12) << 4 |
                       (n.sw->content.b_l & 3) << 8 |
                       (n.se->content.b_l & 3) << 10 |
                       (n.sw->content.b_l & 12) << 10 |
                       (n.se->content.b_l & 12) << 12;
        r = mkblock_leaf(result[unpack_x]);
    } else {
        // Half-sized subblocks of x on the north, south, east, west, and
        // center:
        block *n, *s, *w, *e, *c;
        node no = x->content.b_n;
        assert(no.nw->tag == NODE_B);
        node tmp; // Not an actual node; just a convient way to store four
                  // blocks associated with corners.
        // Recall mkblock_node(nw, ne, sw, se)
        // This part is tedious and error-prone
        n = mkblock_node(no.nw->content.b_n.ne,
                         no.ne->content.b_n.nw,
                         no.nw->content.b_n.se,
                         no.ne->content.b_n.sw);
        s = mkblock_node(no.sw->content.b_n.ne,
                         no.se->content.b_n.nw,
                         no.sw->content.b_n.se,
                         no.se->content.b_n.sw);
        w = mkblock_node(no.nw->content.b_n.sw,
                         no.nw->content.b_n.se,
                         no.sw->content.b_n.nw,
                         no.sw->content.b_n.ne);
        e = mkblock_node(no.ne->content.b_n.sw,
                         no.ne->content.b_n.se,
                         no.se->content.b_n.nw,
                         no.se->content.b_n.ne);
        c = mkblock_node(no.nw->content.b_n.se,
                         no.ne->content.b_n.sw,
                         no.sw->content.b_n.ne,
                         no.se->content.b_n.nw);
        n = evolve(n);
        s = evolve(s);
        w = evolve(w);
        e = evolve(e);
        c = evolve(c);
        tmp.nw = evolve(mkblock_node(evolve(no.nw), n, w, c));
        tmp.ne = evolve(mkblock_node(n, evolve(no.ne), c, e));
        tmp.sw = evolve(mkblock_node(w, c, evolve(no.sw), s));
        tmp.se = evolve(mkblock_node(c, e, s, evolve(no.se)));
        r = mkblock_node(tmp.nw, tmp.ne, tmp.sw, tmp.se);
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
                n.nw = write_bit(n.nw, y, x, bit);
            } else {
                n.ne = write_bit(n.ne, y, x-size/2, bit);
            }
        } else {
            if (x < size/2) {
                n.sw = write_bit(n.sw, y-size/2, x, bit);
            } else {
                n.se = write_bit(n.se, y-size/2, x-size/2, bit);
            }
        }
        if (n.nw == NULL || n.ne == NULL || n.sw == NULL || n.se == NULL) {
            return NULL;
        } else {
            return mkblock_node(n.nw, n.ne, n.sw, n.se);
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
            print_line(b->content.b_n.nw, y, f);
            print_line(b->content.b_n.ne, y, f);
        } else {
            print_line(b->content.b_n.sw, y - size/2, f);
            print_line(b->content.b_n.se, y - size/2, f);
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

main() {
    init_hashtable();
    init_hash_cache();
    init_result();
    int i;
    for (i=0; i<256; i++) {
//        printf("x[%d]: %lu\ny[%d]: %lu\n", i, xmul_cache[i], i, ymul_cache[i]);
    }
    block *b;
    b = read_life_105(stdin);
    if (b == NULL) {
        fprintf(stderr, "Badly formatted input\n");
        exit(1);
    }
    //display(evolve(b), stdout);
    display(b, stdout);
    printf("%lu\n", b->hash);
    //display(b, stdout);
    exit(0);
}
