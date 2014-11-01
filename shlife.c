#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// I'm planning on using bignums for coordinates to allow sizes >2^64. Currently
// this import doesn't have any use except for an almost-pointless use of
// mp_bitcnt_t in the definition of 'block'.
#include <gmp.h>

//// DEFINITIONS OF TYPES

struct block_struct;
typedef struct block_struct block;

// Leaf node for 4x4^H^H^H2x2 block (let's worry about efficiency later).
// Don't change this macro; this value is hardcoded in other places.
#define LEAFSIZE 2
typedef uint8_t leaf;

// A node block recursively containing 4 smaller blocks.
typedef struct {
    block *nw, *ne, *sw, *se;
} node;

// General block. May be a leaf or node.
// To be defined: subblocks of bigger blocks.
struct block_struct {
    enum b_tag {
        EMPTY, // For uninitialized block entries in the hash table.
        LEAF_B,
        NODE_B
        // , CONTAIN_B
    } tag;
    union {
        leaf b_l;
        node b_n;
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

// Note: hashprime must be less than 2^31 for code to work. See comments in
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

int
init_hash_cache() {
    // Values xmul^i*ymul^j for 0 =< i,j < LEAFSIZE
    unsigned long point_value[LEAFSIZE * LEAFSIZE];
    int i, j, yhash, xyhash;
    yhash = 1;
    for (i=0; i<LEAFSIZE; i++) {
        xyhash = yhash;
        for (j=0; j<LEAFSIZE; j++) {
            printf("%d: %lu\n", i*LEAFSIZE + j, xyhash);
            point_value[i*LEAFSIZE + j] = xyhash;
            xyhash = xyhash * xmul % hashprime;
        }
        yhash = yhash * ymul % hashprime;
    }

    int p, tmp;
    unsigned long hash;
    for (p=0; p < 1 << (LEAFSIZE*LEAFSIZE); p++) {
        hash = 0;
        tmp = p;
        for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
            if (tmp & 1) {
                hash = (hash + 2*point_value[i]) % hashprime;
            } else {
                hash = (hash + 1*point_value[i]) % hashprime;
            }
            tmp = tmp >> 1;
        }
        leaf_hash_cache[p] = hash;
    }

    hash = point_value[LEAFSIZE-1] * xmul % hashprime;
    printf("x: %lu\n", hash);
    for (i=0; i < 256; i++) {
        xmul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }

    hash = point_value[LEAFSIZE*(LEAFSIZE-1)] * ymul % hashprime;
    printf("y: %lu\n", hash);
    for (i=0; i<256; i++) {
        ymul_cache[i] = hash;
        hash = hash * hash % hashprime;
    }
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
    uint64_t xmul_d, ymul_d, xymul_d, nwh, neh, swh, seh;
    nwh = nw->hash; neh = ne->hash; swh = sw->hash; seh = se->hash;
    xmul_d = xmul_cache[d];
    ymul_d = ymul_cache[d];
    xymul_d = xmul_d * ymul_d % hashprime;
    // Since hashprime is less than 2^31, all the terms in the sum are less than
    // 2^62, guaranteeing that there is no overflow.
    hash = (nwh + xmul_d*neh + ymul_d*swh + xymul_d*seh) % hashprime;
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
