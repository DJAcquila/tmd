package org.hyperledger.fabric.trustmd.data;

public enum TrustDecision {
    Positive(1),
    Negative(0);

    final private Integer decision;

    TrustDecision(final Integer decision) {
        this.decision = decision;
    }

    public Integer getDecision() {
        return decision;
    }
}