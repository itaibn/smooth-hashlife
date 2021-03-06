#include "shlife.h"
#include "list.c"

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

void mpz_set_size_shift(mpz_t rop, block *b, int i) {
    mpz_set_ui(rop, 1);
    mpz_mul_2exp(rop, rop, LGLENGTH(b)+i);
}

// Two function maybe useful for mulmod

void my_mpz_set_hash(mpz_t rop, const hash_t h) {
    assert (sizeof (hash_t) <= 2 * sizeof (unsigned long));
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(tmp, h >> 32);
    mpz_mul_2exp(rop, tmp, 32);
    mpz_set_ui(tmp, h & (((hash_t) 1>>32)-1));
    mpz_ior(rop, rop, tmp);
    mpz_clear(tmp);
}

hash_t my_mpz_get_hash(mpz_t op) {
    assert (sizeof (hash_t) <= 2 * sizeof (unsigned long));
    assert (mpz_sgn(op) >= 0);
    hash_t low, high;
    mpz_t tmp;
    mpz_init(tmp);
    // mpz_get_ui is documented to return least significant bits of a unsigned
    // value to big to fit.
    low = mpz_get_ui(op) & ((1>>32)-1);
    mpz_tdiv_q_2exp(tmp, op, 32);
    high = mpz_get_ui(tmp);
#ifdef DEBUG
    mpz_tdiv_q_2exp(tmp, tmp, 32);
    assert(mpz_sgn(tmp) == 0);
#endif
    mpz_clear(tmp);
    return (high << 32) | low;
} 

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

// A method for modular multiplication to increase the possible modulo size.
// With this fully incorporated hashprime may be anything less than 2^45.
hash_t mulmod(hash_t a, hash_t b, hash_t m) {
/*
    assert (m < (hash_t) 1 << 45);
    hash_t alow, blow, ahigh, bhigh, res;
    alow = a & ((1<<30)-1);
    blow = b & ((1<<30)-1);
    ahigh = a >> 30;
    bhigh = b >> 30;
    res = 0;
    res += (alow * blow) % m;
res += (((ahigh * blow) % m) << 30) % m;
res += (((alow * bhigh) % m) << 30) % m;
    res += (((ahigh * bhigh) % m) * (((hash_t) 1 << 60) % m)) % m;
    return (res % m);
*/
    mpz_t az, bz, mz, resz;
    hash_t res;
    mpz_inits(az, bz, mz, resz, NULL);
    my_mpz_set_hash(az, a);
    my_mpz_set_hash(bz, b);
    my_mpz_set_hash(mz, m);
    mpz_mul(resz, az, bz);
    mpz_tdiv_r(resz, resz, mz);
    res = my_mpz_get_hash(resz);
    mpz_clears(az, bz, mz, resz, NULL);
    return res;
}

int
init_hash_cache() {
    int i, j, x, y;
    int p, tmp;
    hash_t yhash, xyhash, hash;

    mpz_init_set_ui(hashprime_mpz, hashprime);

    yhash = 1;
    for (y=0; y<LEAFSIZE; y++) {
        xyhash = yhash;
        for (x=0; x<LEAFSIZE; x++) {
            point_hash_value[y*LEAFSIZE + x] = xyhash;
            xyhash = mulmod(xyhash, xmul, hashprime);
        }
        yhash = mulmod(yhash, ymul, hashprime);
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
    hash_t *table, point_value;
    for (y=1; y<=LEAFSIZE; y++) {
    for (x=1; x<=LEAFSIZE; x++) {
        size = 1 << (x*y);
        table = (hash_t *) malloc(size * sizeof(hash_t));
        for (p=0; p < size; p++) {
            hash = 0;
            for (j=0; j<y; j++) {
            for (i=0; i<x; i++) {
                point_value = point_hash_value[LEAFSIZE*j + i];
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

    hash = mulmod(point_hash_value[LEAFSIZE-1], xmul, hashprime);
    for (i=0; i < 256; i++) {
        xmul_cache[i] = hash;
        hash = mulmod(hash, hash, hashprime);
    }

    hash = mulmod(point_hash_value[LEAFSIZE*(LEAFSIZE-1)], ymul, hashprime);
    for (i=0; i<256; i++) {
        ymul_cache[i] = hash;
        hash = mulmod(hash, hash, hashprime);
    }

    return 0;
}

// Note: d is the depth of the subnodes, not the combined node.
hash_t
hash_node(hash_t hnw, hash_t hne, hash_t hsw, hash_t
        hse, depth_t d) {
    if (d >= 256) {
        fprintf(stderr, "This implementation currently does not supported sizes"
            "larger than 2^257\n");
        exit(1);
        //return NULL;
    }

    uint64_t xmul_d, ymul_d, xymul_d;
    hash_t hne_adj, hsw_adj, hse_adj;
    xmul_d = xmul_cache[d];
    hne_adj = mulmod(hne, xmul_d, hashprime);
    ymul_d = ymul_cache[d];
    hsw_adj = mulmod(hsw, ymul_d, hashprime);
    xymul_d = mulmod(xmul_d, ymul_d, hashprime);
    hse_adj = mulmod(hse, xymul_d, hashprime);
    // OLD DOC:
    // Since hashprime is less than 2^31, all the terms in the sum are less than
    // 2^62, guaranteeing that there is no overflow.
    //
    // UPDATE: 
    // With the current code hashprime must be less than 2^62 to ensure no
    // overflow.
    return (hnw + hne_adj + hsw_adj + hse_adj) % hashprime;
}

// Here comes the most complicated function among the hashing utilities.
// Compute the hash of a rectangular subblock with northwest corner (x0, y0) and
// southeast corner (x1, y1). The rectangle is truncated if either (x0, y0) or
// (x1, y1) extend past the ends of base. If the parameter 'adjust' is set
// nonzero, the hash is calculated as if (x0, y0) is the origin.
hash_t
hash_rectangle(block *base, const mpz_t ix0, const mpz_t ix1, const mpz_t iy0,
    const mpz_t iy1, int adjust) {

    hash_t hash;
    mpz_t tmp, blocksize, zero, x0, x1, y0, y1;

    mpz_init_set_ui(zero, 0);
    mpz_init_set_ui(blocksize, LEAFSIZE);
    mpz_mul_2exp(blocksize, blocksize, base->depth);
    mpz_inits(tmp, x0, x1, y0, y1, NULL);
    my_mpz_max(x0, ix0, zero);
    my_mpz_max(y0, iy0, zero);
    my_mpz_min(x1, ix1, blocksize);
    my_mpz_min(y1, iy1, blocksize);
    
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
        hash_t x0l, x1l, y0l, y1l;
        size_t size;
        hash_t *table;
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
        hash = mulmod(point_hash_value[x0l + LEAFSIZE*y0l], table[rect],
            hashprime);
        goto end;
    } else if (base->tag == CONTAIN_B) {
        subblock s = base->content.b_c;
        block *super;
        mpz_t superx0, superx1, supery0, supery1, x_adj, y_adj;
        uint64_t xy_adj;

        super = s.superblock;
        mpz_inits(superx0, superx1, supery0, supery1, x_adj, y_adj, NULL);
        mpz_add(superx0, x0, s.x);
        mpz_add(superx1, x1, s.x);
        mpz_add(supery0, y0, s.y);
        mpz_add(supery1, y1, s.y);

        hash = hash_rectangle(super, superx0, superx1, supery0, supery1, 0);

        mpz_set_ui(x_adj, xmul);
        mpz_neg(tmp, x0);
        mpz_powm(x_adj, x_adj, tmp, hashprime_mpz);
        mpz_init_set_ui(y_adj, ymul);
        mpz_neg(tmp, y0);
        mpz_powm(y_adj, y_adj, tmp, hashprime_mpz);
        xy_adj = mulmod(mpz_get_ui(x_adj), mpz_get_ui(y_adj), hashprime);
        hash = mulmod(xy_adj, hash, hashprime);

        mpz_clears(superx0, superx1, supery0, supery1, NULL);
        goto end;
    }

    assert(base->tag == NODE_B);

    node n = base->content.b_n;
    hash_t hnw, hne, hsw, hse; 
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

end:
    if (adjust) {
        mpz_t x_adj, y_adj;
        uint64_t xy_adj;
        mpz_init_set_ui(x_adj, xmul);
        mpz_neg(tmp, x0);
        mpz_powm(x_adj, x_adj, tmp, hashprime_mpz);
        mpz_init_set_ui(y_adj, ymul);
        mpz_neg(tmp, y0);
        mpz_powm(y_adj, y_adj, tmp, hashprime_mpz);
        xy_adj = mulmod(mpz_get_ui(x_adj), mpz_get_ui(y_adj), hashprime);
        hash = mulmod(xy_adj, hash, hashprime);
        mpz_clears(x_adj, y_adj, NULL);
    }
    mpz_clears(tmp, blocksize, zero, x0, x1, y0, y1, NULL);
    return hash;
}

//// BASIC BLOCK CREATION FUNCTIONS

// Allocate a block with a given hash from the hash table
block * new_block(hash_t hash);

block *
mkblock_leaf(leaf l) {
    hash_t hash = leaf_hash_cache[l];
    block *b = new_block(hash);
    if (b->tag != EMPTY) {
        assert(b->tag == LEAF_B);
        return b;
    }
    b->tag = LEAF_B;
    b->content.b_l = l;
    b->depth = 0;
    return b;
}

// Combine four blocks into a node block
block *
mkblock_node(block *nw, block *ne, block *sw, block *se) {
    hash_t hash;
    block *b;
    depth_t d;

    if (nw == NULL || ne == NULL || sw == NULL || se == NULL) {
        assert(0);
        return NULL;
    }
    d = nw->depth;
    if (ne->depth != d || sw->depth != d || se->depth != d) {
        assert(0);
        return NULL;
    }

    node n = {{{nw, ne}, {sw, se}}};
    hash = hash_node(nw->hash, ne->hash, sw->hash, se->hash, d);
    b = new_block(hash);
    if (b->tag == NODE_B) {
        assert(b->depth == d+1);
        return b;
    }
    if (b->tag == CONTAIN_B) {
        assert(b->depth == d+1);
        if (b->content.b_c.as_node != NULL) {
            return b->content.b_c.as_node;
        } else {
            b->content.b_c.as_node = malloc(sizeof(block));
            b = b->content.b_c.as_node;
            b->hash = hash;
            b->res = NULL;
            b->nfocus = -1;
            b->foci = NULL;
        }
    }
    b->tag = NODE_B;
    b->content.b_n = n;
    b->depth = d+1;
    return b;
}

block *block_index(block *b, int y, int x);

// TODO: Seperate this into two different functions, one for diff=1, one for
// diff>1 (possibly on diff=2 necessary). Thus far this procedure is a mixup
// between block creation and more complicated block processing.
block *
mkblock_contain(block *superblock, mpz_t x, mpz_t y, depth_t diff) {
    assert(superblock);
    assert(mpz_sgn(x) >= 0);
    assert(mpz_sgn(y) >= 0);
    assert(diff >= 0);
    assert(diff <= superblock->depth);

    depth_t d;
    d = superblock->depth;

    if (diff == 0) {
        assert(mpz_sgn(x) == 0 && mpz_sgn(y) == 0);
        return superblock;
    }

    // Verify (x,y) not too far south or east.
    mpz_t rjust, supersize;
    mpz_inits(rjust, supersize, NULL);

    mpz_init_set_ui(supersize, 1);
    mpz_mul_2exp(supersize, supersize, LGLENGTH(superblock));

    mpz_set_ui(rjust, 1);
    mpz_mul_2exp(rjust, rjust, LGLENGTH(superblock) - diff);
    mpz_add(rjust, rjust, x);
    assert (mpz_cmp(rjust, supersize) <= 0);

    mpz_set_ui(rjust, 1);
    mpz_mul_2exp(rjust, rjust, LGLENGTH(superblock) - diff);
    mpz_add(rjust, rjust, y);
    assert (mpz_cmp(rjust, supersize) <= 0);

    mpz_clears(rjust, supersize, NULL);

    if (superblock->tag == CONTAIN_B) {
        mpz_t nx, ny;
        block *res;
        mpz_init(nx); mpz_init(ny);
        mpz_add(nx, x, superblock->content.b_c.x);
        mpz_add(ny, y, superblock->content.b_c.y);
        assert(superblock->content.b_c.superblock);
        res = mkblock_contain(superblock->content.b_c.superblock, nx, ny,
            diff+1);
        mpz_clear(nx); mpz_clear(ny);
        assert(res->depth > 0 || res->tag == LEAF_B);
        return res;
    }
    
    if (diff > 1) {
        mpz_t tmp, nx, ny;
        block *subsuperblock, *res;
        unsigned long x_approx, y_approx;
        depth_t logsize = LGLENGTH(superblock);

        mpz_init(tmp);
        mpz_init(nx); mpz_init(ny);
        mpz_tdiv_q_2exp(tmp, x, logsize-2);
        x_approx = mpz_get_ui(tmp);
        if (x_approx == 3) {
            x_approx = 2;
            mpz_set_ui(tmp, 2);
        }
        mpz_mul_2exp(tmp, tmp, logsize-2);
        mpz_sub(nx, x, tmp);
        mpz_tdiv_q_2exp(tmp, y, logsize-2);
        y_approx = mpz_get_ui(tmp);
        if (y_approx == 3) {
            y_approx = 2;
            mpz_set_ui(tmp, 2);
        }
        mpz_mul_2exp(tmp, tmp, logsize-2);
        mpz_sub(ny, y, tmp);
        mpz_clear(tmp);

        assert(x_approx < 3 && y_approx < 3);
        subsuperblock = block_index(superblock, (int) y_approx, (int) x_approx);
        
        res = mkblock_contain(subsuperblock, nx, ny, diff-1);
        mpz_clear(nx); mpz_clear(ny);
        assert(res->depth > 0 || res->tag == LEAF_B);
        return res;
    }

    assert(diff == 1);
    assert(superblock->tag == NODE_B);
    
    if (d == 1) {
        unsigned long xi, yi;
        xi = mpz_get_ui(x);
        yi = mpz_get_ui(y);
        assert(xi <= 2 && yi <= 2);
        return block_index(superblock, (int) xi, (int) yi);
    }

    hash_t hash;
    mpz_t size, x_east, y_south;
    block *b;

    assert(d >= 0);
    mpz_init_set_ui(size, LEAFSIZE);
    mpz_mul_2exp(size, size, d-diff);
    mpz_init(x_east);
    mpz_add(x_east, x, size);
    mpz_init(y_south);
    mpz_add(y_south, y, size);

    mpz_mul_2exp(size, size, diff);
    assert(mpz_sgn(x) >= 0 && mpz_cmp(x_east, size) <= 0 &&
        mpz_sgn(y) >= 0 && mpz_cmp(y_south, size) <= 0
        && superblock);

    hash = hash_rectangle(superblock, x, x_east, y, y_south, 1);
    mpz_clears(size, x_east, y_south, NULL);

    b = new_block(hash);
    if (b->tag != EMPTY) {return b;}
    b->tag = CONTAIN_B;
    b->depth = d - diff;
    b->content.b_c.superblock = superblock;
    b->content.b_c.as_node = NULL;
    mpz_init_set(b->content.b_c.x, x);
    mpz_init_set(b->content.b_c.y, y);
    // ... b->content.b_c.t ...
    assert(b->depth > 0 || b->tag == LEAF_B);
    return b;
}

subblock
make_subblock_struct(block *superblock, mpz_t x, mpz_t y, depth_t diff) {
    assert(diff >= 1);
    assert(superblock && superblock->tag == NODE_B);

    if (diff == 1) {
        subblock res;
        mpz_init_set(res.y, y);
        mpz_init_set(res.x, x);
        res.superblock = superblock;
        return res;
    }

    mpz_t tmp, nx, ny;
    block *subsuperblock;
    unsigned long x_approx, y_approx;
    depth_t logsize = LGLENGTH(superblock);

    mpz_init(tmp);
    mpz_init(nx); mpz_init(ny);
    mpz_tdiv_q_2exp(tmp, x, logsize-2);
    x_approx = mpz_get_ui(tmp);
    if (x_approx == 3) {
        x_approx = 2;
        mpz_set_ui(tmp, 2);
    }
    mpz_mul_2exp(tmp, tmp, logsize-2);
    mpz_sub(nx, x, tmp);
    mpz_tdiv_q_2exp(tmp, y, logsize-2);
    y_approx = mpz_get_ui(tmp);
    if (y_approx == 3) {
        y_approx = 2;
        mpz_set_ui(tmp, 2);
    }
    mpz_mul_2exp(tmp, tmp, logsize-2);
    mpz_sub(ny, y, tmp);
    mpz_clear(tmp);

    assert(x_approx < 3 && y_approx < 3);
    subsuperblock = block_index(superblock, (int) y_approx, (int) x_approx);
    
    subblock res;
    res = make_subblock_struct(subsuperblock, nx, ny, diff-1);
    mpz_clear(nx); mpz_clear(ny);
    return res;
}

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
block_index(block *b, int iy, int jx) {
    assert(0 <= iy && iy <= 2 && 0 <= jx && jx <= 2);

    if (b == NULL) {return NULL;}
    if (b->tag == LEAF_B) {return NULL;}

    if (b->tag == CONTAIN_B) {
        // Temporarily ad hoc tracing: execeution does not get here.
        assert(0);
        mpz_t x, y, halfblock, shift;
        block *res;

        mpz_inits(x, y, halfblock, shift, NULL);
        mpz_set_ui(halfblock, 1);
        mpz_mul_2exp(halfblock, halfblock, LGLENGTH(b)-2);
        mpz_mul_ui(shift, halfblock, jx);
        mpz_add(x, b->content.b_c.x, shift);
        mpz_mul_ui(shift, halfblock, iy);
        mpz_add(y, b->content.b_c.y, shift);
        res = mkblock_contain(b->content.b_c.superblock, x, y, 2);
        mpz_clears(x, y, halfblock, shift, NULL);
        
        return res;
    } else if (b->tag != NODE_B) {
        fprintf(stderr, "Invalid tag %d for block at %p (hash %lu)\n", b->tag,
            b, b->hash);
        exit(1);
    }

    assert(b->tag == NODE_B);
    node b_no = b->content.b_n;

    if (((iy | jx) & 1) == 0) {
        // This isn't just for optimization; block_index(b, i, j) where either i
        // or j is odd recursively calls the case where both i and j are even.
        return CORNER(b_no, iy/2, jx/2);
    } else if (b->depth == 1) {
    // The next two cases have similar structures. The block b is determined to
    // be a node. The subblock is computed by computing its (p, q)'th corner as
    // (p, q) varies from (0,0) to (1,1). i0 is the higher order bit of p+i and
    // i1 is the lower order bit, and similarly (j0,j1) are the bits of q+j.
    //
    // The two cases are for when the b has depth one and the result is a leaf
    // and when b has depth greater than one and the result is a node. Perhaps
    // the proper way to do this is to have a 'subleaf' block which consists of
    // a quarter of a leaf as well as a general method for combining blocks into
    // nodes as well as subleafs into leafs.

        leaf res, bit;
        int p, i0, i1, q, j0, j1;
        res = 0;
        for (p=0; p<2; p++) {
        for (q=0; q<2; q++) {
            i0 = ((p + iy) & 2) >> 1;
            i1 = (p + iy) & 1;
            j0 = ((q + jx) & 2) >> 1;
            j1 = (q + jx) & 1;
            assert(CORNER(b_no, i0, j0)->tag == LEAF_B);
            bit = 1 & (CORNER(b_no, i0, j0)->content.b_l >> (LEAFSIZE*i1 + j1));
            res |= bit << (LEAFSIZE*p + q);
        }}
        return mkblock_leaf(res);
    } else {
        node n;
        block *tmpb;
        int p, i0, i1, q, j0, j1;
        for (p=0; p<2; p++) {
        for (q=0; q<2; q++) {
            i0 = ((p + iy) & 2) >> 1;
            i1 = (p + iy) & 1;
            j0 = ((q + jx) & 2) >> 1;
            j1 = (q + jx) & 1;
            tmpb = CORNER(b_no, i0, j0);
            assert(tmpb->depth == b->depth - 1);
            assert(((2*i1|2*j1)&1) == 0);
            CORNER(n, p, q) = block_index(tmpb, 2*i1, 2*j1);
        }}
        return mkblock_node(NW(n), NE(n), SW(n), SE(n));
    }
}

//// CA COMPUTATION PROPER

int
copy_inner_pattern(struct inner_pattern *a, struct inner_pattern *b) {
    a->depth_diff = b->depth_diff;
    a->pattern = b->pattern;
    mpz_init_set(a->y, b->y);
    mpz_init_set(a->x, b->x);
    return 1;
}

// Unsure if necessary. Currently unused. [Now used for evolve_with_hint]
block *fill_inner_pattern(block *base, struct inner_pattern *in) {  
    mpz_t north, west, adj;
    mpz_inits(north, west, adj, NULL);
    mpz_set_size_shift(adj, base, -in->depth_diff-1);
    mpz_sub(north, in->y, adj);
    mpz_sub(west, in->x, adj);
    in->pattern = mkblock_contain(base, west, north, in->depth_diff);
    mpz_clears(north, west, adj, NULL);
    return in->pattern;
}

int
is_focal(block *base, struct inner_pattern *compare, int index) {
    struct inner_pattern test = compare[index];
    mpz_t lmargin, rmargin;
    mpz_inits(lmargin, rmargin, NULL);
    mpz_set_size_shift(lmargin, base, -2);
    if ((mpz_cmp(test.x, lmargin) < 0) || (mpz_cmp(test.y, lmargin) < 0)) {
        return 0;
    }
    mpz_set_size_shift(rmargin, base, 0);
    mpz_sub(rmargin, rmargin, lmargin);
    if ((mpz_cmp(test.x, rmargin) > 0) || (mpz_cmp(test.y, rmargin) > 0)) {
        return 0;
    }
    mpz_clears(lmargin, rmargin, NULL);

    struct inner_pattern *comparand = compare;
    mpz_t dist, maxdist;
    mpz_inits(dist, maxdist, NULL);
    mpz_set_size_shift(maxdist, base, -3);

    int i=0;
    while (comparand->depth_diff >= 0) {
        i++;
        if (comparand == &compare[index]) {
            goto continue0;
        }
        if (comparand->pattern->hash < test.pattern->hash) {
            goto continue0;
        }
        mpz_sub(dist, comparand->x, test.x);
        mpz_abs(dist, dist);
        if (mpz_cmp(dist, maxdist) > 0) {
            goto continue0;
        }
        mpz_sub(dist, comparand->y, test.y);
        mpz_abs(dist, dist);
        if (mpz_cmp(dist, maxdist) > 0) {
            goto continue0;
        }
        mpz_clears(dist, maxdist, NULL);
        return 0;

        // Slight modification of 'continue' since that is a keyword.
        continue0:
        comparand++;
    }

    mpz_clears(dist, maxdist, NULL);
    return 1;
}

int
add_foci_node(block *b) {
    assert(b && b->depth >= 3);
    assert(b->nfocus < 0);

    mpz_t lmargin, rmargin, qsize, shiftx, shifty, tmpx, tmpy;
    mpz_inits(lmargin, rmargin, qsize, shiftx, shifty, tmpx, tmpy, NULL);
    mpz_set_size_shift(qsize, b, -2);
    mpz_set_size_shift(lmargin, b, -3);
    mpz_set_size_shift(rmargin, b, 0);
    mpz_sub(rmargin, rmargin, lmargin);
    struct inner_pattern candidates[999], foci[999];
    foci[0].depth_diff = 0;
    // Avoid a stack overflow
    //struct inner_pattern *candidates, *foci;
    //candidates = malloc(999*sizeof(struct inner_pattern));
    //foci = malloc(999*sizeof(struct inner_pattern));
    int iy, jx, k, count, altcount;
    block *index;

    count = 0;
    for (iy=0; iy<3; iy++) {
    for (jx=0; jx<3; jx++) {
        (mpz_mul_ui(shifty, qsize, iy));
        (mpz_mul_ui(shiftx, qsize, jx));
        assert(b && b->depth);
        index = block_index(b, iy, jx);

        if (add_foci(index) < 0) {
            // Exit instead of returning in failure to avoid memory hole with
            // MPZ types
            fprintf(stderr, "add_foci() failed\n");
            exit(EXIT_FAILURE);
            return -1;
        }

        for (k=0; k<index->nfocus; k++) {
            if (   (iy == 0 && (mpz_cmp (index->foci[k].y, lmargin) < 0))
                || (iy == 2 && (mpz_cmp (index->foci[k].y, rmargin) > 0))
                || (jx == 0 && (mpz_cmp (index->foci[k].x, lmargin) < 0))
                || (jx == 2 && (mpz_cmp (index->foci[k].x, rmargin) > 0))) {
                continue;
            }

            assert(count < 999);
            mpz_init(candidates[count].x);
            mpz_init(candidates[count].y);
            mpz_add(candidates[count].y, index->foci[k].y, shifty);
            mpz_add(candidates[count].x, index->foci[k].x, shiftx);
            candidates[count].depth_diff = 2;
            // Ugly: lmargin is exactly the shift necessary to go from the
            // center to the NW corner of the focal candidate.
            mpz_sub(tmpy, candidates[count].y, lmargin);
            mpz_sub(tmpx, candidates[count].x, lmargin);
            candidates[count].pattern = mkblock_contain(b, tmpx, tmpy, 2);
            if (count == 0) {
            }
            count++;
        }
    }}

    mpz_clears(lmargin, rmargin, qsize, shiftx, shifty, tmpx, tmpy, NULL);

    assert(count < 999);
    candidates[count].depth_diff = -1;

    altcount = 0;

//    assert(mpz_sgn(candidates[0].x) > 0);
    for (k=0; k<count; k++) {
        if (is_focal(b, candidates, k)) {
            assert(altcount < 999 && k < 999);
            copy_inner_pattern(&foci[altcount], &candidates[k]);
            altcount++;
        }
    }
    for (k=0; k<count; k++) {
        mpz_clears(candidates[k].x, candidates[k].y, NULL);
    }

    b->nfocus = altcount;
    b->foci = malloc(b->nfocus * sizeof(struct inner_pattern));
    for (k=0; k<altcount; k++) {
        assert(k < b->nfocus && k < 999);
        copy_inner_pattern(&b->foci[k], &foci[k]);
        mpz_clears(foci[k].x, foci[k].y, NULL);
    }

    //free(candidates); free(foci);
    return b->nfocus;
}

// This code was written from scratch and not tested so surely fails badly. It
// also could use some better organizing, I think.
int
add_foci(block *b) {
    assert(b);
    if (b->nfocus >= 0) {
        return b->nfocus;
    }


    block *tmp=NULL;
    // TODO: Determine necessary length for *tmpfoci{0,1}
    //struct inner_pattern tmpfoci0[999], tmpfoci1[999];
    int i, j, k, count;
    count = 0;

    if (b->depth < 3) {
        if (b->depth < 2) {
            b->foci = NULL;
            return (b->nfocus = 0);
        }

        struct inner_pattern *tmpfoci0 = malloc(25 * sizeof(struct
            inner_pattern));
        mpz_t iz, jz;
        mpz_inits(iz, jz, NULL);
        for (j=2; j<7; j++) {
        for (i=2; i<7; i++) {
            mpz_set_ui(jz, j-1);
            mpz_set_ui(iz, i-1);
            tmpfoci0[count].pattern = mkblock_contain(b, iz, jz, 2);
            tmpfoci0[count].depth_diff = 2;
            mpz_inits(tmpfoci0[count].y, tmpfoci0[count].x, NULL);
            mpz_add_ui(tmpfoci0[count].y, jz, 1);
            mpz_add_ui(tmpfoci0[count].x, iz, 1);
            count++;
        }}
        mpz_clears(iz, jz, NULL);

        b->nfocus = count;
        b->foci = (struct inner_pattern *) malloc(count * sizeof(struct
            inner_pattern));
        for (i=0; i<count; i++) {
            b->foci[i] = tmpfoci0[i];
        }
        free(tmpfoci0);
    } else if (b->tag == CONTAIN_B) {
        struct inner_pattern tmpfoci0[999];
        mpz_t mnorth, msouth, mwest, meast, lhs, hsize;
        mpz_inits(mnorth, msouth, meast, mwest, lhs, hsize, NULL);

        mpz_set_ui(mnorth, 1);
        mpz_mul_2exp(hsize, mnorth, LGLENGTH(b));
        mpz_mul_2exp(mnorth, mnorth, LGLENGTH(b)-2); // Am I sure it should be 2
        mpz_set(mwest, mnorth);                      // here?
        mpz_sub(msouth, hsize, mnorth);
        mpz_sub(meast, hsize, mnorth);               // PS. This is a horrible
        mpz_add(mnorth, mnorth, b->content.b_c.y);   // way to format comments.
        mpz_add(msouth, msouth, b->content.b_c.y);
        mpz_add(mwest, mwest, b->content.b_c.x);
        mpz_add(meast, meast, b->content.b_c.x);
        mpz_tdiv_q_2exp(hsize, hsize, 1);

        for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            tmp = block_index(b->content.b_c.superblock, i, j);
            if (add_foci(tmp) < 0) {return -1; /*MPZ NOT DEALLOCATED!!!*/}
            for (k=0; k<tmp->nfocus; k++) {
                int cond;
                /*
                fx, fy = tmp->foci[k].{x,y} + {i,j}
                x, y = b->content.b_c.{x,y}
                cond = fx + j*length(b) >= x + length(b)/4
                    && fx + j*length(b) <= x + length(b)*3/4
                    && fy + i*length(b) >= y + length(b)/4
                    && fy + i*length(b) <= y + length(b)*3/4
                */
                mpz_mul_ui(lhs, hsize, i);
                mpz_add(lhs, lhs, tmp->foci[k].y);
                cond = (mpz_cmp (lhs, mnorth) >= 0) && (mpz_cmp (lhs, msouth)
                    <= 0);
                mpz_mul_ui(lhs, hsize, j);
                mpz_add(lhs, lhs, tmp->foci[k].x);
                cond = cond && (mpz_cmp (lhs, mwest) >= 0) && (mpz_cmp (lhs,
                    meast) <= 0);
                if (cond) {
//                    foci[count++] = tmp->foci[k];
                    tmpfoci0[count].pattern = tmp->foci[k].pattern;
                    tmpfoci0[count].depth_diff = tmp->foci[k].depth_diff;
                    mpz_init(tmpfoci0[count].y);
                    mpz_mul_ui(lhs, hsize, j);
                    mpz_add(tmpfoci0[count].y, tmp->foci[k].y, lhs);
                    mpz_sub(tmpfoci0[count].y, tmpfoci0[count].y, b->content.b_c.y);
                    mpz_init(tmpfoci0[count].x);
                    mpz_mul_ui(lhs, hsize, i);
                    mpz_add(tmpfoci0[count].x, tmp->foci[k].x, lhs);
                    mpz_sub(tmpfoci0[count].x, tmpfoci0[count].x,
                        b->content.b_c.x);
                    count++;
                }
            }
        }}
        mpz_clears(mnorth, msouth, meast, mwest, lhs, hsize, NULL);

        b->nfocus = count;
        b->foci = (struct inner_pattern *) malloc(count * sizeof(struct
            inner_pattern));
        for (i=0; i<count; i++) {
            b->foci[i] = tmpfoci0[i];
        }
    } else if (b->tag == NODE_B) {
        return add_foci_node(b);
    }
    return b->nfocus;
}

int
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

    return 0;
}

block *evolve(block *x);

block *
half_evolve(block *b, int i, int j) {
    assert(b);
    assert(0 <= i && i < 2 && 0 <= j && j < 2);

    node n;
    block *tmp;
    int x,y;

    for (x=0; x<2; x++) {
    for (y=0; y<2; y++) {
        tmp = evolve(block_index(b, i+x, j+y));
        assert(tmp != NULL);
        CORNER(n, x, y) = tmp;
    }}

    return mkblock_node(NW(n), NE(n), SW(n), SE(n));
}

block *evolve_subblock(subblock sb) {
    /*
    mpz_t tmp, newx, newy;
    block *new_sblock;
    int x_approx, y_approx;
    depth_t lglen = LGLENGTH(sb.superblock)-1;

    mpz_inits(tmp, newx, newy);
    mpz_tdiv_q_2exp(tmp, sb.x, lglen);
    x_approx = (int) mpz_get_ui(tmp);
    mpz_mul_2exp(tmp, tmp, lglen);
    mpz_tdiv_q_2exp(tmp, sb.y, lglen);
    y_approx = (int) mpz_get_ui(tmp);
    mpz_clears(tmp, newx, newy);

    assert((x_approx | y_approx) < 2);
    */
    block *res, *new_sblock;
    mpz_t newy, newx, threshold;
    int ix, jy;
    mpz_inits(newy, newx, threshold, NULL);
    mpz_set_size_shift(threshold, sb.superblock, -2);

    if (mpz_cmp(sb.y, threshold) > 0) {
        mpz_sub(newy, sb.y, threshold);
        jy = 1;
    } else {
        mpz_set(newy, sb.y);
        jy = 0;
    }
    if (mpz_cmp(sb.x, threshold) > 0) {
        mpz_sub(newx, sb.x, threshold);
        ix = 1;
    } else {
        mpz_set(newx, sb.x);
        ix = 0;
    }

    assert(sb.superblock->tag == NODE_B);
    new_sblock = half_evolve(sb.superblock, ix, jy);
    res = mkblock_contain(new_sblock, newx, newy, 1);
    mpz_clears(newy, newx, threshold, NULL);
    return res;
}

block *evolve_with_hint(block *b, struct inner_pattern hint);

block *
evolve(block *x) {
    static int nevolve = 0;
    static int cache = 0;
    static int success;
    if (x->depth > 3) {
    }
    nevolve++;
    block *r;
    if (x == NULL) {
        fprintf(stderr, "NULL input to evolve()\n");
    }
    if (x->res) {
        cache++;
        return x->res;
    }
    if (x->tag == CONTAIN_B) {
        r = evolve_subblock(x->content.b_c);
        goto end;
// HOW THE HECK DID THIS EVER WORK IN THE FIRST PLACE?!?!?!
/*
        mpz_t tmp;
        block *new_sblock;
        unsigned long x_approx, y_approx;

        mpz_init(tmp);
        mpz_tdiv_q_2exp(tmp, x->content.b_c.x, LGLENGTH(x));
        x_approx = mpz_get_ui(tmp);
        mpz_tdiv_q_2exp(tmp, x->content.b_c.y, LGLENGTH(x));
        y_approx = mpz_get_ui(tmp);
        mpz_clear(tmp);

        assert((x_approx | y_approx) < 2);

        assert(x->content.b_c.superblock->tag == NODE_B);
        new_sblock = half_evolve(x->content.b_c.superblock, (int) x_approx,
            (int) y_approx);
        r = mkblock_contain(new_sblock, x->content.b_c.x, x->content.b_c.y, 1);
        goto end;
*/
    }
    if (x->tag != NODE_B) {
        fprintf(stderr, "Only nodes have evolve() implemented\n");
        return NULL;
    }
    if (x->depth == 1) {
        // Subblocks are all leafs
        node n = x->content.b_n;
        assert(NW(n)->tag == LEAF_B);
        int unpack_x = (NW(n)->content.b_l & 3) |
                       ((NE(n)->content.b_l & 3) << 2) |
                       ((NW(n)->content.b_l & 12) << 2) |
                       ((NE(n)->content.b_l & 12) << 4) |
                       ((SW(n)->content.b_l & 3) << 8) |
                       ((SE(n)->content.b_l & 3) << 10) |
                       ((SW(n)->content.b_l & 12) << 10) |
                       ((SE(n)->content.b_l & 12) << 12);
        r = mkblock_leaf(result[unpack_x]);
    } else {
        if (add_foci(x) > 0) {
            return evolve_with_hint(x, x->foci[0]);
        }

        int i, j;
        node n;
        assert(x->tag == NODE_B);
        for (i=0; i<2; i++) {
        for (j=0; j<2; j++) {
            CORNER(n, i, j) = evolve(half_evolve(x, i, j));
        }}
        r = mkblock_node(NW(n), NE(n), SW(n), SE(n));
    }
    end:
    if (x->depth > 9) {
        TRACE("evolve %d %d %d %d\n", x->depth, nevolve, cache, success);
    }
    return (x->res = r);
}

// Rough draft
int
intersect(block *block, struct inner_pattern other) {
    mpz_t block_size, shift, north, south, west, east;
    mpz_inits(block_size, shift, north, south, west, east, NULL);
    mpz_set_size_shift(block_size, block, 0);
    mpz_set_size_shift(shift, block, -other.depth_diff-1);
    mpz_sub(north, other.y, shift);
    mpz_add(south, other.y, shift);
    mpz_sub(west, other.x, shift);
    mpz_add(east, other.x, shift);
    /*
    if (mpz_cmp(other.y, block_size) >= 0
     || mpz_cmp(other.x, block_size) >= 0) {
        res = 0;
        goto end;
    }
    if (mpz_cmp(
    res = 1;
    end:
    */
    int res;
    res = (mpz_sgn(south) > 0)
        && (mpz_sgn(east) > 0)
        && (mpz_cmp(north, block_size) < 0)
        && (mpz_cmp(west, block_size) < 0);
    mpz_clears(block_size, shift, north, south, west, east, NULL);
    return res;
}

block *half_evolve_with_hint(block *b, int iy, int jx, struct inner_pattern
    hint);

block *
evolve_with_hint(block *b, struct inner_pattern hint) {
    if (b->depth < 2 || b->tag == CONTAIN_B) {return evolve(b);}

    if (hint.depth_diff < 0 && (mpz_sgn(hint.x) <= 0) && (mpz_sgn(hint.y) <= 0))
    {
        mpz_t low_lim, b_size;
        mpz_inits(low_lim, b_size, NULL);
        mpz_set_size_shift(low_lim, b, -hint.depth_diff);
        mpz_set_size_shift(b_size, b, 0);
        mpz_sub(low_lim, low_lim, b_size);
        mpz_neg(low_lim, low_lim);
        if ((mpz_cmp(low_lim, hint.x) <= 0) && (mpz_cmp(low_lim, hint.y) <= 0))
        {
            mpz_clears(low_lim, b_size, NULL);
            // b is completely contained in hint.pattern
            TRACE("Found subblock %d\n", (int) b->depth);
            subblock sb = make_subblock_struct(hint.pattern, hint.x, hint.y,
                -hint.depth_diff);
            return b->res = evolve_subblock(sb);
        }
    }

    if (!intersect(b, hint)) {
        return evolve(b);
    }

    assert(b->tag == NODE_B);

    int iy, jx;
    node n;
    block *half_evolve;
    //struct inner_pattern new_hint;
    //mpz_t shift, shift_mul;
    
    /*
    mpz_inits(new_hint.x, new_hint.y, shift, shift_mul, NULL);
    new_hint.depth_diff = hint.depth_diff - 1;
    new_hint.pattern = hint.pattern;
    mpz_set_size_shift(shift_mul, b, -3);
    */
    for (iy=0; iy<2; iy++) {
    for (jx=0; jx<2; jx++) {
        half_evolve = half_evolve_with_hint(b, iy, jx, hint);
        assert(half_evolve->depth == b->depth - 1);
        /*
        mpz_mul_ui(shift, shift_mul, 1+2*iy);
        mpz_sub(new_hint.y, hint.y, shift);
        mpz_mul_ui(shift, shift_mul, 1+2*jx);
        mpz_sub(new_hint.x, hint.x, shift);
        */
        CORNER(n, iy, jx) = evolve(half_evolve);
    }}
    //mpz_clears(new_hint.x, new_hint.y, shift, shift_mul, NULL);

    b->res = mkblock_node(NW(n), NE(n), SW(n), SE(n));
    assert(b->res->depth == b->depth-1);
    return b->res;
}

block *
half_evolve_with_hint(block *b, int iy, int jx, struct inner_pattern hint) {
    assert(b);
    assert(0 <= iy && iy < 2 && 0 <= jx && jx < 2);

    node n;
    block *tmp;
    int x,y;

    for (x=0; x<2; x++) {
    for (y=0; y<2; y++) {
        // ??? This shouldn't work !!!
        tmp = evolve_with_hint(block_index(b, iy+x, jx+y), hint);
        assert(tmp != NULL && tmp->depth == b->depth-2);
        CORNER(n, x, y) = tmp;
    }}

    block *res = mkblock_node(NW(n), NE(n), SW(n), SE(n));
    assert(res->depth == b->depth-1);
    return res;
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

int display_raw(block *b, FILE *f);

block *
blank_block(depth_t lglength) {
    if (lglength < LGLEAFSIZE) {
        fprintf(stderr, "too small argument to blank_block()\n");
    }

    block *blank;
    blank = mkblock_leaf(0);
    while (LGLENGTH(blank) < lglength) {
        blank = mkblock_node(blank, blank, blank, blank);
    }
    return blank;
}

// To aid a lazy implementation, a maximum length for lines
#define MAXMCLINES 1000000
block *
read_mc(FILE *f) {

    char c;

    // Read first line: "[M2] ..."
    while ((c = getc(f)) != '\n') {
        if (c == EOF) {
            return NULL;
        }
    } 
    // Read second line: "#..."
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
    blocktable = malloc(MAXMCLINES * sizeof(block *));
    blocktable[0] = NULL;
    blocktable[1] = NULL;
    //if (DEBUG) display_raw(blank_eight);

    //int inw, ine, isw, ise, depth;
    int depth;
    int i,j;
    int index;
    node n;
    index = 1;
    while(index < MAXMCLINES) {
        /*
        if (DEBUG) {
            //TRACE("previous patterns (index=%d):\n", index);
            for (i=0; i<index; i++) {
                //TRACE("(%d)", i);
                display_raw(blocktable[i], stderr);
            }
            //TRACE("\n");
        }
        */
        block *current;
        c = getc(f);
        if (c == '$' || c == '*' || c == '.') {
            int x=0, y=0, bit;
            current = blank_eight;
            while (c != '\n') {
                switch(c) {
                    case '.':
                        bit = 0;
                        goto dot_or_star;
                    case '*':
                        bit = 1;
                    dot_or_star:
                        current = write_bit(current, y, x, bit);
                        assert(x++ < 8);
                        break;
                    case '$':
                        x = 0;
                        assert(y++ < 8);
                        break;
                    default:
                        fprintf(stderr, "Invalid formatting in Macrocell "
                            "file\n");
                        return NULL; // WARNING: Doesn't deallocate resources.
                }
                c = getc(f);
            }
        } else if ('0' <= c && c <= '9') {
            int x[5];
            for (i=0; i<5; i++) {
                x[i] = 0;
                do {
                    assert(x[i] < INT_MAX / 10);
                    x[i] = x[i]*10 + (c - '0');
                    c = getc(f);
                } while ('0' <= c && c <= '9');
                if (i<4) {
                    assert(c == ' ');
                    c = getc(f);
                    assert('0' <= c && c <= '9');
                }
            }
            if (c == '\r') {c = getc(f);}
            assert(c == '\n');
            assert(x[1] < index && x[2] < index && x[3] < index && x[4] <
                index);

            depth = x[0];
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                if (x[2*i+j+1] == 0) {
                    CORNER(n, i, j) = blank_block(depth-1);
                } else if (x[2*i+j+1] < index) {
                    CORNER(n, i, j) = blocktable[x[2*i+j+1]];
                } else {
                    fprintf(stderr, "Bad formatting\n");
                    return NULL; // WARNING: Doesn't deallocate resources.
                }
            }}
            current = mkblock_node(NW(n), NE(n), SW(n), SE(n));
            assert(current && LGLENGTH(current) == depth);
        } else if (c == EOF) {
            break;
        } else {
            fprintf(stderr, "Nonsensical line in .mc");
            exit(EXIT_FAILURE);
        }

        blocktable[index] = current;
        index++;
    }

    TRACE("done\n");
    return blocktable[index-1];
}

// Useful for debugging.
int
display_raw(block *b, FILE *f) {
    if (!b) {
        fprintf(f, "NULL");
        return 1;
    }

    node n;
    block *btmp;
    switch (b->tag) {
        case LEAF_B:
            fprintf(f, "%x", b->content.b_l);
            break;
        case CONTAIN_B:
            gmp_fprintf(f, "c[%Zd %Zd]", b->content.b_c.x, b->content.b_c.y);
            display_raw(b->content.b_c.superblock, f);
            break;
        case NODE_B:
            fprintf(f, "n(");
            int i, j;
            n = b->content.b_n;
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                btmp = CORNER(n, i, j);
                display_raw(btmp, f);
            }}
            fprintf(f, ")");
            break;
        default:
            fprintf(stderr, "Invalid tag for block input for display()\n");
            return 1;
    }
    return 0;
}

// Used in read_life_105. Returns NULL when (x,y) is out of range.
block *
write_bit(block *b, unsigned long y, unsigned long x, char bit) {
    unsigned long size = 2;
    size = 2 << b->depth;

    if (x >= size) {
        return NULL;
    }
    if (y >= size) {
        return NULL;
    }

    if (b->tag == LEAF_B) {
        leaf mask = 1 << (2*y+x);
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

// Display the line from (x0, y) to (x1, y) in Life 1.05 format onto file f.
// Ends with a newline if newline!=0.
int
display_line(block *b, mpz_t y, mpz_t x0, mpz_t x1, int newline, FILE *f) {
    assert(mpz_cmp(x0, x1) <= 0);
    mpz_t size, halfsize, tmpx0, tmpx1, tmpy;

    mpz_init(size); mpz_init(halfsize); mpz_init(tmpx0); mpz_init(tmpx1);
    mpz_init(tmpy);
    mpz_set_ui(size, LEAFSIZE);
    mpz_mul_2exp(size, size, b->depth);
    mpz_tdiv_q_2exp(halfsize, size, 1);

    if (mpz_sgn(y) < 0 || mpz_cmp(y, size) >= 0 || mpz_sgn(x1) <= 0 ||
        mpz_cmp(x0, size) >= 0) {
        goto end;
    }

    unsigned long yl, x0l, x1l;
    int i, j, bit;
    leaf bt;
    node nb;
    subblock sb;
    switch (b->tag) {
        case LEAF_B:
            yl = mpz_get_ui(y);
            if (mpz_sgn(x0) < 0) {
                x0l = 0;
            } else {
                x0l = mpz_get_ui(x0);
            }
            if (mpz_cmp(x1, size) > 0) {
                x1l = LEAFSIZE;
            } else {
                x1l = mpz_get_ui(x1);
            }
            bt = b->content.b_l;

            assert(x0l < x1l);
            for (i=x0l; i<x1l; i++) {
                bit = 1 & (bt >> (yl*LEAFSIZE + i));
                fputc(bit ? '*' : '.', f);
            }
            break;

        case NODE_B:
            nb = b->content.b_n;
            for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                mpz_mul_si(tmpy, halfsize, -i);
                mpz_add(tmpy, tmpy, y);
                mpz_mul_si(tmpx0, halfsize, -j);
                mpz_add(tmpx0, tmpx0, x0);
                mpz_mul_si(tmpx1, halfsize, -j);
                mpz_add(tmpx1, tmpx1, x1);

                display_line(CORNER(nb, i, j), tmpy, tmpx0, tmpx1, 0, f);
            }}
            break;

        case CONTAIN_B:
            sb = b->content.b_c;

            mpz_set_ui(tmpy, 0);
            my_mpz_max(tmpx0, x0, tmpy);
            my_mpz_min(tmpx1, x1, size);

            mpz_add(tmpy, y, sb.y);
            mpz_add(tmpx0, tmpx0, sb.x);
            mpz_add(tmpx1, tmpx1, sb.x);

            display_line(sb.superblock, tmpy, tmpx0, tmpx1, 0, f);
            break;

        default:
            fprintf(stderr, "Invalid block type\n");
            return 1;
    }

    end:
    if (newline) {
        fputc('\n', f);
    }

    mpz_clear(size); mpz_clear(halfsize); mpz_clear(tmpx0); mpz_clear(tmpx1);
    mpz_clear(tmpy);
    return 0;
}

// Displays a block in Life 1.05 format.
int
display(block *b, FILE *f) {
    int ret=0;
    int size;
    size = LEAFSIZE << b->depth;
    if (size == 0) {
        fprintf(stderr, "Can't display block: Too big\n");
        return 1;
    }

    //TRACE("d size %d", size);
    //display_raw(b, stdout);
    //TRACE("\n");

    fputs("#Life 1.05\n", f);

    int i;
    mpz_t y, x0, x1;
    mpz_init_set_ui(y, 0);
    mpz_init_set_ui(x0, 0);
    mpz_init_set_ui(x1, size);
    for (i=0; i<size; i++) {
        //TRACE("d l %Zd ", y);
        if (display_line(b, y, x0, x1, 1, f)) {
            fprintf(stderr, "Error display line %d of block\n", i);
            ret = 1;
            goto end;
        }
        mpz_add_ui(y, y, 1);
    }
end:
    mpz_clears(x0, x1, y, NULL);
    return ret;
}

//// HASH TABLE-RELATED FUNCTIONS
// Note: Currently this sucks, there is no hash table resizing nor garbage
// collection.

/*
unsigned int ht_size = 1000000;
block **hashtable;
*/

int
init_hashtable() {
    hashtable = (block **) malloc(ht_size * sizeof(block *));
    memset(hashtable, 0, ht_size * sizeof(block *));
    return 0;
}

block *
new_block(hash_t hash) {
    hash_t ii = 0;
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
            b->nfocus = -1;
            b->foci = NULL;
#if DEBUG
            b->index = i;
#endif
//            b->refcount = 0;
            // Insert command here of the form:
            //  num_entries++;
            hashtable[i] = b;
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

int
initialize() {
    return (init_hashtable() || init_hash_cache() || init_result());
}

void
show_point_hash_cache() {
    int i;
    printf("Point hash cache:\n");
    for (i=0; i<LEAFSIZE*LEAFSIZE; i++) {
        printf("%d: %llu\n", i, point_hash_value[i]);
    }
}

void
show_leaf_hash_cache() {
    int i;
    printf("Leaf hash cache:\n");
    for (i=0; i<1<<(LEAFSIZE*LEAFSIZE); i++) {
        printf("%x: %llu\n", i, (unsigned long long) leaf_hash_cache[i]);
    }
}

int
main(int argc, char** argv) {
    init_hashtable();
    init_hash_cache();
    init_result();

    if (argc != 3) {
        fprintf(stderr, "shlife pattern.mc out.lif\n");
        exit(1);
    }

    FILE *f;
    block *b;
    f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(1);
    }
    b = read_mc(f);
    //b = read_life_105(f);
    if (b == NULL) {
        fprintf(stderr, "Badly formatted input\n");
        exit(1);
    }
    close(f);
    assert(sizeof(unsigned long long) == sizeof(hash_t));
    show_point_hash_cache();
    show_leaf_hash_cache();
    printf("hash: %llu\n", (unsigned long long) b->hash);
    printf("mulmod test: %llu\n", (unsigned long long) mulmod((hash_t)1<<35,
        (hash_t)1<<35, hashprime));

    f = fopen(argv[2], "w");
    if (b->depth > 0) {
        display(evolve(b), f);
    }
    close(f);
}

/*
int
main(int argc, char **argv) {
    //TRACE("EMPTY %d LEAF_B %d NODE_B %d CONTAIN_B %d\n", EMPTY, LEAF_B, NODE_B,
    //    CONTAIN_B);

    init_hashtable();
    init_hash_cache();
    init_result();

    if (argc != 7) {
        fprintf(stderr, "shlife pattern outfile x0 x1 y0 y1\n");
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
    //b = read_mc(f);
    if (b == NULL) {
        fprintf(stderr, "Badly formatted input\n");
        exit(1);
    } else if (b->depth == 0) {
        fprintf(stderr, "Input pattern must be larger than leaf\n");
        // Not exiting intentionally.
    }
    fclose(f);

    mpz_t x[4];
    int i, j, err=0;
    for (i=0; i<4; i++) {
        err |= mpz_init_set_str(x[i], argv[i+3], 10);
    }
    if (err) {
        fprintf(stderr, "Arguments x0, x1, y0, y1 must be numbers\n");
        exit(1);
    }

    f = fopen(argv[2], "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening output file\n");
        exit(1);
    }

    block *b1 = block_index(b, 0, 1);
    mpz_t three, five;
    mpz_init(three); mpz_init(five);
    mpz_set_ui(three, 3);
    mpz_set_ui(five, 5);
    //TRACE("TEST mkc\n");
    //mkblock_contain(b1, three, three, 2, 0);
    
    int n;
    TRACE("FOCALS d %d\n", b->depth);
    printf("Focals: %d\n", n = add_foci(b));
    display(b, stdout);
    for (i=0; i<n; i++) {
        struct inner_pattern in = b->foci[i];
        gmp_printf("(%Zd, %Zd), ", in.x, in.y);
    }
    printf("\n");
    printf("Focals of variant:\n");
    block *d = blank_block(LGLENGTH(b));
    //node no;
    //NW(no) = b; NE(no) = d; SW(no) = d; SE(no) = d;
    d = mkblock_node(b, d, d, d);
    d = mkblock_contain(d, x[0], x[1], 1);
    printf("Focals: %d\n", n = add_foci(d));
    display(d, stdout);
    for (i=0; i<n; i++) {
        struct inner_pattern in = d->foci[i];
        gmp_printf("(%Zd, %Zd), ", in.x, in.y);
    }
    printf("\n");


    printf("evolving\n");
    b = evolve(b);
    printf("done\n");
    if (b == NULL) {
        fprintf(stderr, "Error evolving life pattern\n");
        exit(1);
    }
    display(b, f);
    fclose(f);

    exit(0);
}
*/
