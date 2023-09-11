package org.trustmd.model;

import com.owlike.genson.Genson;
import org.eclipse.microprofile.openapi.annotations.media.Schema;

@Schema(name="Trust")
public class Trust {

    private final static Genson genson = new Genson();
    private String mecHostId;
    private String clusterHeadId;
    private String evaluatedNodeId;
    private String trustDecision;
    private String degreeOfBelief;

    public Trust(final String clusterHeadId,
                 final String evaluatedNodeId, final String trustDecision,
                 final String degreeOfBelief, final String mecHostId) {
        this.clusterHeadId = clusterHeadId;
        this.degreeOfBelief = degreeOfBelief;
        this.trustDecision = trustDecision;
        this.mecHostId = mecHostId;
        this.evaluatedNodeId = evaluatedNodeId;
    }

    public Trust() {

    }

    /**
     * FUNC : Serialize the Class and its Data.
     *
     * @return class : String | Serilaized Class of the Tag.
     */
    public String toJSONString() {
        return genson.serialize(this);
    }

    /**
     * FUNC : Deserialize a string and create an instance of a NEW Tag
     *
     * @param json : String | JSON string to be Deserialized.
     * @return asset : Tag | Newly created object of the Tag.
     */
    public static Trust fromJSONString(String json) {
        return genson.deserialize(json, Trust.class);
    }

    public String getTrustDecision() {
        return trustDecision;
    }

    public void setTrustDecision(final String trustDecision) {
        this.trustDecision = trustDecision;
    }

    public String getDegreeOfBelief() {
        return degreeOfBelief;
    }

    public void setDegreeOfBelief(final String degreeOfBelief) {
        this.degreeOfBelief = degreeOfBelief;
    }

    public String getEvaluatedNodeId() {
        return evaluatedNodeId;
    }

    public void setEvaluatedNodeId(final String evaluatedNodeId) {
        this.evaluatedNodeId = evaluatedNodeId;
    }

    public String getClusterHeadId() {
        return clusterHeadId;
    }

    public void setClusterHeadId(final String clusterHeadId) {
        this.clusterHeadId = clusterHeadId;
    }

    public String getMecHostId() {
        return mecHostId;
    }

    public void setMecHostId(final String mecHostId) {
        this.mecHostId = mecHostId;
    }
}
