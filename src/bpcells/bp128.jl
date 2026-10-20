# Scalar decoder for the BP128 bit-packing used by BPCells "packed" matrices.
#
# Layout (from BPCells `simd/bp128/bp128-helper-inl.h`): 128 unsigned 32-bit values are packed
# with `bits` bits each into `4 * bits` 32-bit words, as four interleaved lanes. Lane `l`
# (0..3) holds values `4i + l` for `i = 0..31`, packed LSB-first into the word sequence
# `words[4k + l]`, `k = 0..bits-1`. `bits == 0` means all zero; `bits == 32` is a plain copy.
# The D1 variant stores deltas: value[j] = start + sum(delta[0..j]); D1Z zigzag-encodes the
# deltas first (so a drop at a column boundary stays small). BPCells packed matrices use D1Z
# for indices and FOR (frame 1) for integer values.

"""
    bp128_unpack!(out, words, woff, bits)

Unpack one 128-value chunk into `out` (length ≥ 128) from `words[woff+1 : woff+4*bits]`.
"""
function bp128_unpack!(out::AbstractVector{UInt32}, words::AbstractVector{UInt32}, woff::Int, bits::Int)
    if bits == 0
        fill!(view(out, 1:128), UInt32(0))
        return out
    elseif bits == 32
        copyto!(out, 1, words, woff + 1, 128)
        return out
    end
    mask = (UInt32(1) << bits) - UInt32(1)
    @inbounds for lane in 0:3
        k = 0                       # next word (within this lane) to load
        inreg = UInt32(0)
        for i in 0:31
            shift = (i * bits) & 31
            if shift == 0
                inreg = words[woff + 4k + lane + 1]; k += 1
                outreg = inreg
            else
                outreg = inreg >> shift
            end
            if shift > 32 - bits
                inreg = words[woff + 4k + lane + 1]; k += 1
                outreg |= inreg << (32 - shift)
            end
            out[4i + lane + 1] = outreg & mask
        end
    end
    return out
end

"""Unpack a D1 (delta-coded) chunk: cumulative sum of the deltas starting from `start`."""
function bp128_unpack_d1!(out::AbstractVector{UInt32}, words::AbstractVector{UInt32}, woff::Int, bits::Int, start::UInt32)
    bp128_unpack!(out, words, woff, bits)
    acc = start
    @inbounds for j in 1:128
        acc += out[j]
        out[j] = acc
    end
    return out
end

"""Unpack a FOR chunk with frame of reference `frame` (BPCells uses 1 for integer values)."""
function bp128_unpack_for!(out::AbstractVector{UInt32}, words::AbstractVector{UInt32}, woff::Int, bits::Int, frame::UInt32)
    bp128_unpack!(out, words, woff, bits)
    @inbounds for j in 1:128
        out[j] += frame
    end
    return out
end

"""Unpack a D1Z chunk: zigzag-decoded deltas, cumulative from `start` (BPCells packed matrix indices)."""
function bp128_unpack_d1z!(out::AbstractVector{UInt32}, words::AbstractVector{UInt32}, woff::Int, bits::Int, start::UInt32)
    bp128_unpack!(out, words, woff, bits)
    acc = start
    @inbounds for j in 1:128
        u = out[j]
        delta = (u >> 1) ⊻ (-(u & UInt32(1)))     # zigzag decode, wraps modulo 2^32
        acc += delta
        out[j] = acc
    end
    return out
end
