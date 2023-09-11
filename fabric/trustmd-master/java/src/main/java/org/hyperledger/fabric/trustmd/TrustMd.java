/*
Copyright IBM Corp., DTCC All Rights Reserved.

SPDX-License-Identifier: Apache-2.0
*/
package org.hyperledger.fabric.trustmd;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.google.protobuf.ByteString;
import com.owlike.genson.Genson;
import io.netty.handler.ssl.OpenSsl;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hyperledger.fabric.shim.ChaincodeBase;
import org.hyperledger.fabric.shim.ChaincodeException;
import org.hyperledger.fabric.shim.ChaincodeStub;
import org.hyperledger.fabric.shim.ledger.KeyModification;
import org.hyperledger.fabric.shim.ledger.KeyValue;
import org.hyperledger.fabric.shim.ledger.QueryResultsIterator;
import org.hyperledger.fabric.trustmd.data.Trust;

import static java.nio.charset.StandardCharsets.UTF_8;

public class TrustMd extends ChaincodeBase {

    private static final Log _logger = LogFactory.getLog(TrustMd.class);
    private final static Genson genson = new Genson();

    @Override
    public Response init(ChaincodeStub stub) {
        try {
            _logger.info("Init java TrustMdContract");
            final String func = stub.getFunction();
            if (!func.equals("init")) {
                return newErrorResponse("function other than init is not supported");
            }
            final List<String> args = stub.getParameters();
            if (args.size() != 5) {
                return newErrorResponse("Incorrect number of arguments. Expecting 5");
            }

            final Trust trustAsset = buildTrustAsset(args);

            stub.putStringState(trustAsset.getEvaluatedNodeId(), trustAsset.toJSONString());

            return newSuccessResponse();
        } catch (Exception e) {
            return newErrorResponse(e);
        }
    }

    private Trust buildTrustAsset(final List<String> args) {

        final String clusterHeadId = args.get(0);
        final String evaluatedNodeId = args.get(1);
        final String trustDecision =  args.get(2);
        final String degreeOfBelief = args.get(3);
        final String mecHostId = args.get(4);

        _logger.info("Creating trust asset");
        // Initialize the chaincode
        final Trust asset = new Trust(clusterHeadId,
                evaluatedNodeId, trustDecision,
                degreeOfBelief, mecHostId);

        _logger.info(String.format("clusterHeadId %s; evaluatedNodeId %s, trustDecision %s, degreeOfBelief %s, mecHostId %s",
                clusterHeadId, evaluatedNodeId, trustDecision, degreeOfBelief, mecHostId)
        );

        return asset;
    }

    @Override
    public Response invoke(final ChaincodeStub stub) {
        try {
            _logger.info("Invoke java TrustMdContract chaincode");
            final String func = stub.getFunction();
            final List<String> params = stub.getParameters();
            if (func.equals("createTrust")) {
                return createTrust(stub, params);
            }
            if (func.equals("updateTrust")) {
                return updateTrustDecisionAndValue(stub, params);
            }
            if (func.equals("getTrust")) {
                return getTrust(stub, params);
            }
            if (func.equals("delete")) {
                return delete(stub, params);
            }
            return newErrorResponse("Invalid invoke function name. Expecting one of: [\"createTrust\", \"updateTrust\", \"getTrust\", \"getTrust\"]");
        } catch (Throwable e) {
            return newErrorResponse(e);
        }
    }
    public boolean trustExists(final ChaincodeStub stub, final String evaluatedNodeId) {
        final byte[] buffer = stub.getState(evaluatedNodeId);
        return (buffer != null && buffer.length > 0);
    }

    private Response createTrust(final ChaincodeStub stub, final List<String> args) {
        final Trust trustAsset = buildTrustAsset(args);
        final boolean exists = trustExists(stub, trustAsset.getEvaluatedNodeId());
        if (exists) {
            return newSuccessResponse("The asset " + trustAsset.getEvaluatedNodeId() + " already exists");
        }
        try {
            stub.putStringState(trustAsset.getEvaluatedNodeId(), trustAsset.toJSONString());
            return newSuccessResponse("createTrust finished successfully", ByteString.copyFrom(trustAsset.toJSONString().getBytes(UTF_8)).toByteArray());

        } catch (final Exception e) {
            System.out.println(e.toString());
            throw new ChaincodeException("Error", e.toString());
        }
    }

    public Response queryTrustHistory(final ChaincodeStub stub, final List<String> args) {
        if (args.size() != 1) {
            return newErrorResponse("Incorrect number of arguments. Expecting name of the person to query");
        }

        final String evaluatedNodeId = args.get(0);

        if (evaluatedNodeId == null) {
            return newErrorResponse("No evaluatedNodeId given");
        }

        final ArrayList<String> results = new ArrayList<>();
        try {
            final QueryResultsIterator<KeyModification> history = stub.getHistoryForKey(evaluatedNodeId);

            if (history == null) {
                final String errorMessage = String.format("Trust %s does not exist", evaluatedNodeId);
                System.out.println(errorMessage);
                return newErrorResponse(errorMessage);
            }

            for (final KeyModification keyModification : history) {
                final String iteratorValue = keyModification.getStringValue();
                results.add(iteratorValue);
                System.out.println(iteratorValue);
            }
            history.close();
        } catch (Exception e) {
            results.add(e.getMessage());
            results.add(e.getCause().getMessage());
            results.add(Arrays.toString(e.getStackTrace()));
        }
        return newSuccessResponse();
    }

    public Response getTrust(final ChaincodeStub stub, final List<String> args) {
        final String evaluatedNodeId = args.get(0);

        final String trustState = stub.getStringState(evaluatedNodeId);

        if (trustState.isEmpty()) {
            final String errorMessage = String.format("Trust of node %s does not exist", evaluatedNodeId);
            System.out.println(errorMessage);
            return newSuccessResponse(errorMessage);
        }

        return newSuccessResponse(trustState, ByteString.copyFrom(trustState, UTF_8).toByteArray());
    }

    private Response updateTrustDecisionAndValue(final ChaincodeStub stub, final List<String> args) {
        final String evaluatedNodeId = args.get(1);
        final boolean exists = trustExists(stub, evaluatedNodeId);
        if (!exists) {
            return createTrust(stub, args);
        }
        try {
            final Trust newTrust = buildTrustAsset(args);

            stub.putStringState(evaluatedNodeId, newTrust.toJSONString());
            return newSuccessResponse(newTrust.toJSONString(), ByteString.copyFrom(newTrust.toJSONString(), UTF_8).toByteArray());
        } catch (Exception e) {
            System.out.println(e.toString());
            return newErrorResponse("[Error] " + e.toString());
        }
    }

    public Response getAllTrust(final ChaincodeStub stub) {
        final QueryResultsIterator<KeyValue> queryIt = stub.getStateByRange("", "");
        final ArrayList<String> results = new ArrayList<>();
        if (queryIt == null) {
            return newErrorResponse("Empty ChainCode");
        }
        try {
            for (KeyValue keyValue : queryIt) {
                String key = keyValue.getKey();
                results.add(key);
                System.out.println(key);
            }
            queryIt.close();
        } catch (Exception e) {
            results.add(e.getMessage());
            results.add(e.getCause().getMessage());
            results.add(Arrays.toString(e.getStackTrace()));
        }
        return newSuccessResponse(results.toString(), ByteString.copyFrom(results.toString(), UTF_8).toByteArray());
    }

    private Response delete(ChaincodeStub stub, List<String> args) {
        if (args.size() != 1) {
            return newErrorResponse("Incorrect number of arguments. Expecting 1");
        }
        String key = args.get(0);
        // Delete the key from the state in ledger
        stub.delState(key);
        return newSuccessResponse();
    }

    public static void main(String[] args) {
        System.out.println("OpenSSL avaliable: " + OpenSsl.isAvailable());
        new TrustMd().start(args);
    }

}
