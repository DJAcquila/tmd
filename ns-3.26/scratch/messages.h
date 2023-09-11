#include "ns3/ipv4-global-routing-helper.h"

#define CH_CHANGE 0
#define CH_NEW 1
#define CH_NOT_CHANGE 2
#define ADD_MEMBER 3
#define SEARCH_ERR 4

using namespace ns3;

enum Messages
{
    SEND_TRANSACTION,
    VALIDATE_TRANSACTION,
    VOTE,
    FINISH,
    ACK_VOTE,
    ACK_CH,
    START_VOTE,
    CLOSE,  
    GET_NEW_BLOCK,
    BLOCK,
    BLOCK_VOTE,
	INV,        
	GET_HEADERS, 
	HEADERS,   
	GET_BLOCKS,       
	GET_DATA,   
	NO_MESSAGE, 
	SEND_PUBLIC_KEY,
	RECEIVE_PUBLIC_EDSA_KEY, 
    RECEIVE_PUBLIC_EDH_KEY,
	RECEIVE_MESSAGE  
};

enum BlockBroadcastType
{
	STANDARD             //DEFAULT
};

enum ProtocolType
{
	STANDARD_PROTOCOL,             //DEFAULT
	SENDHEADERS
};

typedef struct {
    int nodeId;
    int manufacturerId;
    Ipv4Address ipv4Address;
    std::string nodePublicKey;
    std::string signature;
} blockDataTuple;

typedef struct {
        double downloadSpeed;
        double uploadSpeed;
} nodeInternetSpeeds;


const char* getMessageName(enum Messages m)
{
  switch (m)
  {
    case SEND_TRANSACTION: return "SEND_TRANSACTION";
    case VALIDATE_TRANSACTION: return "VALIDATE_TRANSACTION";
    case VOTE: return "VOTE";
    case ACK_VOTE: return "ACK_VOTE";
    case ACK_CH: return "ACK_CH";
    case CLOSE: return "CLOSE";
    case GET_NEW_BLOCK: return "GET_NEW_BLOCK";
    case INV: return "INV";
    case GET_HEADERS: return "GET_HEADERS";
    case HEADERS: return "HEADERS";
    case GET_BLOCKS: return "GET_BLOCKS";
    case BLOCK: return "BLOCK";
    case GET_DATA: return "GET_DATA";
    case NO_MESSAGE: return "NO_MESSAGE";
    case SEND_PUBLIC_KEY: return "SEND_PUBLIC_KEY";
    case RECEIVE_PUBLIC_EDSA_KEY: return "RECEIVE_PUBLIC_EDSA_KEY";
    case RECEIVE_PUBLIC_EDH_KEY: return "RECEIVE_PUBLIC_EDH_KEY";
    case RECEIVE_MESSAGE: return "RECEIVE_MESSAGE";
    case BLOCK_VOTE: return "BLOCK_VOTE";
    case START_VOTE: return "START_VOTE";
    case FINISH: return "FINISH";
  }
  return "ERR";
}

const char* getBlockBroadcastType(enum BlockBroadcastType m)
{
  switch (m)
  {
    case STANDARD: return "STANDARD";
    /*case UNSOLICITED: return "UNSOLICITED";
    case RELAY_NETWORK: return "RELAY_NETWORK";
    case UNSOLICITED_RELAY_NETWORK: return "UNSOLICITED_RELAY_NETWORK";*/
  }
  return "ERR";
}

const char* getProtocolType(enum ProtocolType m)
{
  switch (m)
  {
    case STANDARD_PROTOCOL: return "STANDARD_PROTOCOL";
    case SENDHEADERS: return "SENDHEADERS";
  }
  return "ERR";
}


typedef struct {
        int nodeId;
        bool ch;
        int totalBlocks;
        // int attackSuccess;                   //0->fail, 1->success
        long invReceivedBytes;
        long invSentBytes;
        long getHeadersReceivedBytes;
        long getHeadersSentBytes;
        long headersReceivedBytes;
        long headersSentBytes;
        long getDataReceivedBytes;
        long getDataSentBytes;
        long blockReceivedBytes;
        long blockSentBytes;
        long blockTimeouts;
        double avg_consensusTime;
        // long chunkTimeouts;
        // int minedBlocksInMainChain;
} nodeStatistics;